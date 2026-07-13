# SPDX-License-Identifier: GPL-3.0-or-later

# Boundary-condition impedance. Ports BCImpedancecMod.f90 (complex, KRAKENC),
# BCImpedanceMod.f90 (real, KRAKEN) and misc/PekRoot.f90.

"""
Pekeris branch of the square root exposing leaky modes (ports `PekerisRoot`,
misc/PekRoot.f90).
"""
@inline pekeris_root(z::Complex) = real(z) >= 0 ? sqrt(z) : im * sqrt(-z)

"""
    bc_impedance(w, prob, x, bottop, hs, cnt) -> (f, g, ipower)

Impedance functions (f, g) of the boundary for trial eigenvalue `x` (= k²).
`bottop` is `:top` or `:bot`. The complex method ports `BCImpedance`
(BCImpedancecMod.f90); the real method ports BCImpedanceMod.f90 (`ComplexFlag`
= .FALSE. — the only value the reachable code paths use) and increments the
mode counter `cnt` for an elastic bottom halfspace, as the Fortran does.
"""
function bc_impedance(w::MeshWorkspace{T,Complex{T}}, prob::Problem{T},
                      x::Complex{T}, bottop::Symbol, hs::Halfspace{T},
                      cnt) where {T}
    ω² = (2π * prob.freq)^2
    ipower = 0
    yV = zeros(Complex{T}, 5)

    local f::Complex{T}, g::Complex{T}
    if hs.bc === :vacuum
        f = one(Complex{T})
        g = zero(Complex{T})
        yV[1] = f; yV[2] = g
    elseif hs.bc === :rigid
        f = zero(Complex{T})
        g = one(Complex{T})
        yV[1] = f; yV[2] = g
    else                                    # :acoustoelastic
        if real(hs.cs) > 0                  # elastic halfspace
            gammaS2 = x - ω² / hs.cs^2
            gammaP2 = x - ω² / hs.cp^2
            gammaS = pekeris_root(gammaS2)
            gammaP = pekeris_root(gammaP2)
            mu = hs.rho * hs.cs^2
            yV[1] = (gammaS * gammaP - x) / mu
            yV[2] = ((gammaS2 + x)^2 - 4 * gammaS * gammaP * x) * mu
            yV[3] = 2 * gammaS * gammaP - gammaS2 - x
            yV[4] = gammaP * (x - gammaS2)
            yV[5] = gammaS * (gammaS2 - x)
            f = ω² * yV[4]
            g = yV[2]
        else                                # acoustic halfspace
            gammaP = pekeris_root(x - ω² / hs.cp^2)
            f = gammaP
            g = Complex{T}(hs.rho)
        end
    end

    bottop === :top && (g = -g)    # top BC has the sign flipped vs bottom

    # Shoot through elastic layers
    nmedia = length(prob.media)
    if bottop === :top && w.first_acoustic > 1
        for medium in 1:w.first_acoustic-1
            ipower = elastic_dn!(w, x, yV, ipower, medium)
        end
        f = ω² * yV[4]
        g = yV[2]
    elseif bottop === :bot && w.last_acoustic < nmedia
        for medium in nmedia:-1:w.last_acoustic+1
            ipower = elastic_up!(w, x, yV, ipower, medium)
        end
        f = ω² * yV[4]
        g = yV[2]
    end

    f, g, ipower
end

function bc_impedance(w::MeshWorkspace{T,T}, prob::Problem{T}, x::T,
                      bottop::Symbol, hs::Halfspace{T}, cnt) where {T}
    ω² = (2π * prob.freq)^2
    ipower = 0
    yV = zeros(T, 5)

    local f::Complex{T}, g::Complex{T}
    if hs.bc === :vacuum
        f = one(Complex{T}); g = zero(Complex{T})
        yV[1] = real(f); yV[2] = real(g)
    elseif hs.bc === :rigid
        f = zero(Complex{T}); g = one(Complex{T})
        yV[1] = real(f); yV[2] = real(g)
    else                                    # :acoustoelastic
        if real(hs.cs) > 0                  # elastic halfspace (real parts)
            gammaS2 = x - ω² / real(hs.cs)^2
            gammaP2 = x - ω² / real(hs.cp)^2
            gammaS = real(sqrt(complex(gammaS2)))
            gammaP = real(sqrt(complex(gammaP2)))
            mu = hs.rho * real(hs.cs)^2
            yV[1] = (gammaS * gammaP - x) / mu
            yV[2] = (gammaS2 + x)^2 - 4 * gammaS * gammaP * x
            yV[2] *= mu
            yV[3] = 2 * gammaS * gammaP - gammaS2 - x
            yV[4] = gammaP * (x - gammaS2)
            yV[5] = gammaS * (gammaS2 - x)
            f = Complex{T}(ω² * yV[4])
            g = Complex{T}(yV[2])
            real(g) > 0 && (cnt[] += 1)     # mode count contribution
        else                                # acoustic halfspace
            gammaP = sqrt(x - ω² / hs.cp^2)
            # ComplexFlag = .FALSE.: real parts only
            f = Complex{T}(real(gammaP))
            g = Complex{T}(hs.rho)
        end
    end

    bottop === :top && (g = -g)

    nmedia = length(prob.media)
    if bottop === :top && w.first_acoustic > 1
        for medium in 1:w.first_acoustic-1
            ipower = elastic_dn!(w, x, yV, ipower, medium)
        end
        f = Complex{T}(ω² * yV[4])
        g = Complex{T}(yV[2])
    elseif bottop === :bot && w.last_acoustic < nmedia
        for medium in nmedia:-1:w.last_acoustic+1
            ipower = elastic_up!(w, x, yV, ipower, medium)
        end
        f = Complex{T}(ω² * yV[4])
        g = Complex{T}(yV[2])
    end

    f, g, ipower
end

"""
Propagate through an elastic layer, shooting up from below, using the compound
matrix formulation (ports `ElasticUP`, BCImpedance(c)Mod.f90). Mutates `yV`,
returns the updated power-of-ten exponent.
"""
function elastic_up!(w::MeshWorkspace{T,X}, x, yV::Vector, ipower::Int,
                     medium::Int) where {T,X}
    two_x = 2 * x
    two_h = 2 * w.h[medium]
    four_h_x = 4 * w.h[medium] * x
    j = w.loc[medium] + w.n[medium] + 1
    xB3 = x * w.b3[j] - w.rho[j]

    xV = similar(yV); zV = similar(yV)
    # Euler's method for first step
    zV[1] = yV[1] - (w.b1[j] * yV[4] - w.b2[j] * yV[5]) / 2
    zV[2] = yV[2] - (-w.rho[j] * yV[4] - xB3 * yV[5]) / 2
    zV[3] = yV[3] - (two_h * yV[4] + w.b4[j] * yV[5]) / 2
    zV[4] = yV[4] - (xB3 * yV[1] + w.b2[j] * yV[2] - two_x * w.b4[j] * yV[3]) / 2
    zV[5] = yV[5] - (w.rho[j] * yV[1] - w.b1[j] * yV[2] - four_h_x * yV[3]) / 2

    # Modified midpoint method
    for ii in w.n[medium]:-1:1
        j -= 1
        copyto!(xV, yV)
        copyto!(yV, zV)
        xB3 = x * w.b3[j] - w.rho[j]
        zV[1] = xV[1] - (w.b1[j] * yV[4] - w.b2[j] * yV[5])
        zV[2] = xV[2] - (-w.rho[j] * yV[4] - xB3 * yV[5])
        zV[3] = xV[3] - (two_h * yV[4] + w.b4[j] * yV[5])
        zV[4] = xV[4] - (xB3 * yV[1] + w.b2[j] * yV[2] - two_x * w.b4[j] * yV[3])
        zV[5] = xV[5] - (w.rho[j] * yV[1] - w.b1[j] * yV[2] - four_h_x * yV[3])
        if ii != 1     # scale if necessary
            if abs(real(zV[2])) < FLOOR_
                zV .*= ROOF; yV .*= ROOF
                ipower -= IPOW_R
            end
            if abs(real(zV[2])) > ROOF
                zV .*= FLOOR_; yV .*= FLOOR_
                ipower -= IPOW_F
            end
        end
    end

    yV .= (xV .+ 2 .* yV .+ zV) ./ 4   # standard filter at the terminal point
    ipower
end

"Ports `ElasticDN` (BCImpedance(c)Mod.f90) — shooting down from above."
function elastic_dn!(w::MeshWorkspace{T,X}, x, yV::Vector, ipower::Int,
                     medium::Int) where {T,X}
    two_x = 2 * x
    two_h = 2 * w.h[medium]
    four_h_x = 4 * w.h[medium] * x
    j = w.loc[medium] + 1
    xB3 = x * w.b3[j] - w.rho[j]

    xV = similar(yV); zV = similar(yV)
    zV[1] = yV[1] + (w.b1[j] * yV[4] - w.b2[j] * yV[5]) / 2
    zV[2] = yV[2] + (-w.rho[j] * yV[4] - xB3 * yV[5]) / 2
    zV[3] = yV[3] + (two_h * yV[4] + w.b4[j] * yV[5]) / 2
    zV[4] = yV[4] + (xB3 * yV[1] + w.b2[j] * yV[2] - two_x * w.b4[j] * yV[3]) / 2
    zV[5] = yV[5] + (w.rho[j] * yV[1] - w.b1[j] * yV[2] - four_h_x * yV[3]) / 2

    for ii in 1:w.n[medium]
        j += 1
        copyto!(xV, yV)
        copyto!(yV, zV)
        xB3 = x * w.b3[j] - w.rho[j]
        zV[1] = xV[1] + (w.b1[j] * yV[4] - w.b2[j] * yV[5])
        zV[2] = xV[2] + (-w.rho[j] * yV[4] - xB3 * yV[5])
        zV[3] = xV[3] + (two_h * yV[4] + w.b4[j] * yV[5])
        zV[4] = xV[4] + (xB3 * yV[1] + w.b2[j] * yV[2] - two_x * w.b4[j] * yV[3])
        zV[5] = xV[5] + (w.rho[j] * yV[1] - w.b1[j] * yV[2] - four_h_x * yV[3])
        if ii != w.n[medium]
            if abs(real(zV[2])) < FLOOR_
                zV .*= ROOF; yV .*= ROOF
                ipower -= IPOW_R
            end
            if abs(real(zV[2])) > ROOF
                zV .*= FLOOR_; yV .*= FLOOR_
                ipower -= IPOW_F
            end
        end
    end

    yV .= (xV .+ 2 .* yV .+ zV) ./ 4
    ipower
end
