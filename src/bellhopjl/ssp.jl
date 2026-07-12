# SPDX-License-Identifier: GPL-3.0-or-later

# Sound-speed profiles + derivatives. Ports the depth-dependent interpolants of
# sspMod.f90 (cLinear, n2Linear, cCubic, UpdateDepthSegmentT) and the de Boor
# spline coefficients from misc/splinec.f90 (CSPLINE, end conditions
# iBCBeg = iBCEnd = 0, i.e. not-a-knot).

"Result of an SSP evaluation (range-independent SSPs ⇒ cr = crr = crz = 0)."
struct SSPEval{T<:Real}
    c::T
    cimag::T
    cz::T
    czz::T
end

SSPEval(c, cimag, cz, czz) = SSPEval(promote(c, cimag, cz, czz)...)

abstract type AbstractSSP{T<:Real} end

Base.eltype(::AbstractSSP{T}) where {T} = T

"""
    update_seg(z, zq, tz, iseg) -> iseg′

Depth-segment update with the tangent-direction-aware edge-case handling of
`UpdateDepthSegmentT` (sspMod.f90): if the ray takes a small step in direction
`tz` it remains in the same segment.
"""
@inline function update_seg(z::AbstractVector, zq, tz, iseg::Int)
    n = length(z)
    if tz >= 0
        while zq < z[iseg] && iseg > 1
            iseg -= 1
        end
        while zq >= z[iseg + 1] && iseg < n - 1
            iseg += 1
        end
    else
        while zq > z[iseg + 1] && iseg < n - 1
            iseg += 1
        end
        while zq <= z[iseg] && iseg > 1
            iseg -= 1
        end
    end
    iseg
end

########################### C-linear ('C') ###########################

"C-linear SSP. Ports `cLinear` (sspMod.f90); node speeds are complex (CRCI)."
struct CLinearSSP{T<:Real} <: AbstractSSP{T}
    z::Vector{T}
    c::Vector{Complex{T}}
    cz::Vector{Complex{T}}       # per-segment gradient
end

function CLinearSSP(z::AbstractVector, c::AbstractVector)
    T = float(promote_type(eltype(z), real(eltype(c))))
    zv = collect(T, z)
    cv = collect(Complex{T}, c)
    issorted(zv) || throw(ArgumentError("SSP depths must be increasing"))
    length(zv) ≥ 2 || throw(ArgumentError("SSP needs at least 2 points"))
    CLinearSSP{T}(zv, cv, diff(cv) ./ diff(zv))
end

@inline function evaluate(s::CLinearSSP{T}, zq, tz, iseg::Int) where {T}
    i = update_seg(s.z, zq, tz, iseg)
    cc = s.c[i] + (zq - s.z[i]) * s.cz[i]
    SSPEval(real(cc), imag(cc), real(s.cz[i]), zero(real(cc))), i
end

########################### N²-linear ('N') ###########################

"N²-linear SSP. Ports `n2Linear` (sspMod.f90)."
struct N2LinearSSP{T<:Real} <: AbstractSSP{T}
    z::Vector{T}
    n2::Vector{Complex{T}}
    n2z::Vector{Complex{T}}
end

function N2LinearSSP(z::AbstractVector, c::AbstractVector)
    T = float(promote_type(eltype(z), real(eltype(c))))
    zv = collect(T, z)
    cv = collect(Complex{T}, c)
    issorted(zv) || throw(ArgumentError("SSP depths must be increasing"))
    n2 = 1 ./ cv .^ 2
    N2LinearSSP{T}(zv, n2, diff(n2) ./ diff(zv))
end

@inline function evaluate(s::N2LinearSSP{T}, zq, tz, iseg::Int) where {T}
    i = update_seg(s.z, zq, tz, iseg)
    w = (zq - s.z[i]) / (s.z[i+1] - s.z[i])
    cc = 1 / sqrt((1 - w) * s.n2[i] + w * s.n2[i+1])
    c = real(cc)
    cz = -c^3 * real(s.n2z[i]) / 2
    czz = 3 * cz * cz / c
    SSPEval(c, imag(cc), cz, czz), i
end

########################### Cubic spline ('S') ###########################

"""
Cubic-spline SSP. Ports `cCubic` (sspMod.f90) with de Boor `CSPLINE`
coefficients (misc/splinec.f90), end conditions iBCBeg = iBCEnd = 0
(not-a-knot). In segment i: f(h) = C1 + h(C2 + h(C3/2 + h C4/6)),
h = z - z[i].
"""
struct SplineSSP{T<:Real} <: AbstractSSP{T}
    z::Vector{T}
    coef::Matrix{Complex{T}}      # 4 × n
end

"""
    cspline_coefs(z, y) -> 4×n matrix

Ports CSPLINE (misc/splinec.f90) for the iBCBeg = iBCEnd = 0 (not-a-knot)
case actually used by Bellhop (`cCubic`, sspMod.f90). Direct translation.
"""
function cspline_coefs(tau::Vector{T}, y::Vector{Complex{T}}) where {T}
    n = length(tau)
    C = zeros(Complex{T}, 4, n)
    C[1, :] .= y
    L = n - 1
    for m in 2:n
        C[3, m] = tau[m] - tau[m-1]
        C[4, m] = (C[1, m] - C[1, m-1]) / C[3, m]
    end
    # beginning boundary condition (IBCBEG = 0)
    if n > 2
        C[4, 1] = C[3, 3]
        C[3, 1] = C[3, 2] + C[3, 3]
        C[2, 1] = ((C[3, 2] + 2 * C[3, 1]) * C[4, 2] * C[3, 3] +
                   C[3, 2]^2 * C[4, 3]) / C[3, 1]
    else
        C[4, 1] = 1
        C[3, 1] = 1
        C[2, 1] = 2 * C[4, 2]
    end
    # running calculations to n-1 (not executed if n = 2)
    for m in 2:L
        g = -C[3, m+1] / C[4, m-1]
        C[2, m] = g * C[2, m-1] + 3 * (C[3, m] * C[4, m+1] + C[3, m+1] * C[4, m])
        C[4, m] = g * C[3, m-1] + 2 * (C[3, m] + C[3, m+1])
    end
    # ending boundary condition (IBCEND = 0)
    local g::Complex{T} = 0
    if n == 2       # (N==2 .AND. IBCBEG==0)
        C[2, n] = C[4, n]
    elseif n == 3   # (N==3 .AND. IBCBEG==0)
        C[2, n] = 2 * C[4, n]
        C[4, n] = 1
        g = -1 / C[4, n-1]
    else
        g = C[3, n-1] + C[3, n]
        C[2, n] = ((C[3, n] + 2 * g) * C[4, n] * C[3, n-1] +
                   C[3, n]^2 * (C[1, n-1] - C[1, n-2]) / C[3, n-1]) / g
        g = -g / C[4, n-1]
        C[4, n] = C[3, n-1]
    end
    if n > 2
        C[4, n] = g * C[3, n-1] + C[4, n]
        C[2, n] = (g * C[2, n-1] + C[2, n]) / C[4, n]
    end
    # back substitution
    for j in L:-1:1
        C[2, j] = (C[2, j] - C[3, j] * C[2, j+1]) / C[4, j]
    end
    # final calculations
    for i in 2:n
        dtau = C[3, i]
        divdf1 = (C[1, i] - C[1, i-1]) / dtau
        divdf3 = C[2, i-1] + C[2, i] - 2 * divdf1
        C[3, i-1] = 2 * (divdf1 - C[2, i-1] - divdf3) / dtau
        C[4, i-1] = (divdf3 / dtau) * (6 / dtau)
    end
    # curvature at the last node & mean value (unused by evaluation; kept for parity)
    C[3, n] = C[3, L] + (tau[n] - tau[L]) * C[4, L]
    C[4, n] = 0
    C
end

function SplineSSP(z::AbstractVector, c::AbstractVector)
    T = float(promote_type(eltype(z), real(eltype(c))))
    zv = collect(T, z)
    cv = collect(Complex{T}, c)
    issorted(zv) || throw(ArgumentError("SSP depths must be increasing"))
    length(zv) ≥ 2 || throw(ArgumentError("SSP needs at least 2 points"))
    SplineSSP{T}(zv, cspline_coefs(zv, cv))
end

@inline function evaluate(s::SplineSSP{T}, zq, tz, iseg::Int) where {T}
    i = update_seg(s.z, zq, tz, iseg)
    h = zq - s.z[i]
    # SPLINEALL (misc/splinec.f90)
    f   = s.coef[1, i] + h * (s.coef[2, i] + h * (s.coef[3, i] / 2 + h * s.coef[4, i] / 6))
    fx  = s.coef[2, i] + h * (s.coef[3, i] + h * s.coef[4, i] / 2)
    fxx = s.coef[3, i] + h * s.coef[4, i]
    SSPEval(real(f), imag(f), real(fx), real(fxx)), i
end
