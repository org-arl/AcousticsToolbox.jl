# SPDX-License-Identifier: GPL-3.0-or-later

# Finite-difference mesh setup. Ports `Initialize` (kraken.f90 / krakenc.f90)
# and the SSP subtabulation (`cLinear`, `cCubic`, misc/sspMod.f90).

"""
    cspline_coefs(z, y) -> 4×n matrix

Ports CSPLINE (misc/splinec.f90) for the iBCBeg = iBCEnd = 0 case, as used by
`ReadSSP` (sspMod.f90). Same code as the BellhopJL port (spline SSPs share
this routine upstream).
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
    for m in 2:L
        g = -C[3, m+1] / C[4, m-1]
        C[2, m] = g * C[2, m-1] + 3 * (C[3, m] * C[4, m+1] + C[3, m+1] * C[4, m])
        C[4, m] = g * C[3, m-1] + 2 * (C[3, m] + C[3, m+1])
    end
    local g::Complex{T} = 0
    if n == 2
        C[2, n] = C[4, n]
    elseif n == 3
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
    for j in L:-1:1
        C[2, j] = (C[2, j] - C[3, j] * C[2, j+1]) / C[4, j]
    end
    for i in 2:n
        dtau = C[3, i]
        divdf1 = (C[1, i] - C[1, i-1]) / dtau
        divdf3 = C[2, i-1] + C[2, i] - 2 * divdf1
        C[3, i-1] = 2 * (divdf1 - C[2, i-1] - divdf3) / dtau
        C[4, i-1] = (divdf3 / dtau) * (6 / dtau)
    end
    C[3, n] = C[3, L] + (tau[n] - tau[L]) * C[4, L]
    C[4, n] = 0
    C
end

"""
    tabulate_ssp!(cp, cs, rhov, med, n1)

Subtabulate the medium's SSP at the `n1` mesh points (ports `cLinear` /
`cCubic`, sspMod.f90). Density is rho-linear in both cases (the spline path
splines rho too upstream, but the wrapper writes piecewise-constant rho per
medium, so linear is exact for all reachable inputs — see PORTING_NOTES §3).
"""
function tabulate_ssp!(cp::AbstractVector{Complex{T}}, cs, rhov, med::Medium{T},
                       n1::Int) where {T}
    zn = med.zn
    npts = length(zn)
    n = n1 - 1
    h = (zn[npts] - zn[1]) / n
    cpco = med.spline ? cspline_coefs(zn, med.cpn) : Matrix{Complex{T}}(undef, 0, 0)
    csco = med.spline ? cspline_coefs(zn, med.csn) : Matrix{Complex{T}}(undef, 0, 0)
    lay = 1
    for iz in 1:n1
        z = zn[1] + (iz - 1) * h
        iz == n1 && (z = zn[npts])      # make sure no overshoot
        while lay < npts - 1 && z > zn[lay+1]
            lay += 1
        end
        if med.spline
            hs = z - zn[lay]
            # SPLINE (misc/splinec.f90): f = C1 + h(C2 + h(C3/2 + h C4/6))
            cp[iz] = cpco[1, lay] + hs * (cpco[2, lay] + hs * (cpco[3, lay] / 2 +
                     hs * cpco[4, lay] / 6))
            cs[iz] = csco[1, lay] + hs * (csco[2, lay] + hs * (csco[3, lay] / 2 +
                     hs * csco[4, lay] / 6))
        else
            r = (z - zn[lay]) / (zn[lay+1] - zn[lay])
            cp[iz] = (1 - r) * med.cpn[lay] + r * med.cpn[lay+1]
            cs[iz] = (1 - r) * med.csn[lay] + r * med.csn[lay+1]
        end
        r = (z - zn[lay]) / (zn[lay+1] - zn[lay])
        rhov[iz] = (1 - r) * med.rhon[lay] + r * med.rhon[lay+1]
    end
    nothing
end

"""
    initialize!(w, prob, iset) -> nothing

Set up the finite-difference coefficient arrays for mesh `iset`. Ports
`Initialize` (krakenc.f90 for `X = Complex{T}`, kraken.f90 for `X = T`).
The caller must have set `w.n` and `w.h` (main loop of kraken(c).f90).
"""
function initialize!(w::MeshWorkspace{T,X}, prob::Problem{T}) where {T,X}
    ω² = (2π * prob.freq)^2
    nmedia = length(prob.media)
    w.cmin = T(Inf)
    w.first_acoustic = 0
    w.last_acoustic = 0
    w.loc[1] = 0

    npoints = sum(w.n[1:nmedia]) + nmedia
    resize!(w.b1, npoints);  fill!(w.b1, zero(X))
    resize!(w.b1c, npoints); fill!(w.b1c, zero(T))
    resize!(w.b2, npoints);  fill!(w.b2, zero(X))
    resize!(w.b3, npoints);  fill!(w.b3, zero(X))
    resize!(w.b4, npoints);  fill!(w.b4, zero(X))
    resize!(w.rho, npoints); fill!(w.rho, zero(T))
    cp = Vector{Complex{T}}(undef, npoints)
    cs = Vector{Complex{T}}(undef, npoints)

    elastic_flag = false
    for medium in 1:nmedia
        medium > 1 && (w.loc[medium] = w.loc[medium-1] + w.n[medium-1] + 1)
        n1 = w.n[medium] + 1
        ii = w.loc[medium] + 1

        tabulate_ssp!(view(cp, ii:ii+w.n[medium]), view(cs, ii:ii+w.n[medium]),
                      view(w.rho, ii:ii+w.n[medium]), prob.media[medium], n1)

        if all(iszero, view(cs, ii:ii+w.n[medium]))     # acoustic medium
            w.first_acoustic == 0 && (w.first_acoustic = medium)
            w.last_acoustic = medium
            for j in ii:ii+w.n[medium]
                w.cmin = min(w.cmin, real(cp[j]))
                if X <: Complex
                    w.b1[j] = -2 + w.h[medium]^2 * ω² / cp[j]^2
                else
                    # KRAKEN keeps the real part in B1 and the imaginary part
                    # of ω²/c² in B1C for the attenuation perturbation
                    w.b1[j] = -2 + w.h[medium]^2 * real(ω² / cp[j]^2)
                    w.b1c[j] = imag(ω² / cp[j]^2)
                end
            end
        else                                            # elastic medium
            prob.sigma[medium] == 0 && prob.sigma[medium+1] == 0 ||
                error("Rough elastic interfaces are not allowed")
            elastic_flag = true
            two_h = 2 * w.h[medium]
            for j in ii:ii+w.n[medium]
                w.cmin = min(w.cmin, real(cs[j]))
                cp2 = X <: Complex ? cp[j]^2 : real(cp[j]^2)
                cs2 = X <: Complex ? cs[j]^2 : real(cs[j]^2)
                w.b1[j] = two_h / (w.rho[j] * cs2)
                w.b2[j] = two_h / (w.rho[j] * cp2)
                w.b3[j] = 4 * two_h * w.rho[j] * cs2 * (cp2 - cs2) / cp2
                w.b4[j] = two_h * (cp2 - 2 * cs2) / cp2
                w.rho[j] = two_h * ω² * w.rho[j]
            end
        end
    end

    # (cLow, cHigh) adjustment from the halfspace properties
    chigh = prob.chigh
    if prob.hsbot.bc === :acoustoelastic
        if real(prob.hsbot.cs) > 0
            elastic_flag = true
            w.cmin = min(w.cmin, real(prob.hsbot.cs))
            X <: Complex || (chigh = min(chigh, real(prob.hsbot.cs)))
        else
            w.cmin = min(w.cmin, real(prob.hsbot.cp))
        end
    end
    if prob.hstop.bc === :acoustoelastic
        if real(prob.hstop.cs) > 0
            elastic_flag = true
            w.cmin = min(w.cmin, real(prob.hstop.cs))
            X <: Complex || (chigh = min(chigh, real(prob.hstop.cs)))
        else
            w.cmin = min(w.cmin, real(prob.hstop.cp))
        end
    end

    elastic_flag && (w.cmin *= T(0.85))     # reduce cMin for Scholte wave
    # KRAKENC uses cLow = max(cLow, 0.99 cMin); KRAKEN uses max(cLow, cMin)
    w.clow = X <: Complex ? max(prob.clow, T(0.99) * w.cmin) :
                            max(prob.clow, w.cmin)
    w.chigh = chigh
    nothing
end

"""
Auto mesh count for a medium (ports the NG == 0 rule of ReadEnvironmentMod:
20 points per wavelength at the last-read speed of the medium, min 10).
"""
function auto_ng(med::Medium{T}, freq) where {T}
    c = _value(real(med.cpn[end]))
    iszero(med.csn[end]) || (c = _value(real(med.csn[end])))
    deltaz = c / _value(freq) / 20
    max(floor(Int, _value(med.z1 - med.z0) / deltaz), 10)
end
