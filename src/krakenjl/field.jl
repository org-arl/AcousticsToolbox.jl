# SPDX-License-Identifier: GPL-3.0-or-later

# Pressure field from modes. Ports `Evaluate` (KrakenField/EvaluateMod.f90)
# with the point-source, cylindrical-coordinates option ('R...C'/'R...I')
# that the wrapper's .flp file requests, and `Weight`
# (misc/calculateweights.f90) for mode interpolation at arbitrary depths.

"""
    weight(z, ztab) -> (iz, w)

Linear-interpolation indices and weights of `ztab` into the (ascending) grid
`z` (ports `Weight`, misc/calculateweights.f90).
"""
function weight(z::AbstractVector, ztab::AbstractVector)
    n = length(z)
    iz = similar(ztab, Int)
    w = similar(ztab, Float64)
    for (k, zt) in pairs(ztab)
        i = clamp(searchsortedlast(z, zt), 1, n - 1)
        iz[k] = i
        w[k] = (zt - z[i]) / (z[i+1] - z[i])
    end
    iz, w
end

"Interpolate mode shapes onto depths `ztab` (Vector, kraken(c).f90 PhiTab)."
function tabulate_modes(res::ModeResult{T}, ztab::AbstractVector) where {T}
    iz, wts = weight(res.z, ztab)
    M = length(res.k)
    phitab = Matrix{Complex{T}}(undef, length(ztab), M)
    for m in 1:M, (k, i) in pairs(iz)
        phitab[k, m] = res.phi[i, m] + wts[k] * (res.phi[i+1, m] - res.phi[i, m])
    end
    phitab
end

"""
    evaluate_field(k, phiS, phiR, rxr; incoherent=false, threads=1)
        -> Matrix{Complex} (nz × nr)

Mode summation for a point source in cylindrical coordinates (ports
`Evaluate`, KrakenField/EvaluateMod.f90, Option 'R…C'/'R…I', zero array
tilt). `phiS` are the mode values at the source depth, `phiR` the (nz × M)
mode values at the receiver depths, `rxr` the receiver ranges [m].
Range columns are independent, so they are chunked across threads (each
column written by exactly one task — deterministic).
"""
function evaluate_field(k::Vector{Complex{T}}, phiS::Vector{Complex{T}},
                        phiR::Matrix{Complex{T}}, rxr::AbstractVector;
                        incoherent::Bool=false, threads::Int=1) where {T}
    M = length(k)
    nz = size(phiR, 1)
    nr = length(rxr)
    P = zeros(Complex{T}, nz, nr)
    M <= 0 && return P     # if no modes, return vanishing pressure

    factor = im * sqrt(2 * T(π)) * cis(T(π) / 4)
    cnst = [factor * phiS[m] / sqrt(k[m]) for m in 1:M]     # point source
    ik = [-im * k[m] for m in 1:M]                          # e^{i(ωt - kr)} form
    incoherent && (ik = Complex{T}.(real.(ik)))             # incoherent case
    cmat = permutedims(phiR) .* cnst   # (M × nz): const(m)·phi(m,iz), zero tilt

    hank = [Vector{Complex{T}}(undef, M) for _ in 1:max(threads, 1)]
    fill1range! = (ir, hnk) -> begin
        r = rxr[ir]
        @inbounds for m in 1:M
            hnk[m] = exp(ik[m] * r)
        end
        @inbounds for iz in 1:nz
            s = zero(Complex{T})
            if !incoherent           # coherent case
                for m in 1:M
                    s += cmat[m, iz] * hnk[m]
                end
            else                     # incoherent case (as the Fortran: the
                for m in 1:M         # square, not |·|², under the sqrt)
                    s += (cmat[m, iz] * hnk[m])^2
                end
                s = sqrt(s)
            end
            # cylindrical spreading (Option 'R')
            P[iz, ir] = abs(r) > floatmin(T) ? s / sqrt(r) : s
        end
    end

    if threads <= 1 || nr < 2 * threads
        for ir in 1:nr
            fill1range!(ir, hank[1])
        end
    else
        tasks = map(enumerate(_chunks(nr, threads))) do (t, idxs)
            Threads.@spawn begin
                local hnk = hank[t]
                for ir in idxs
                    fill1range!(ir, hnk)
                end
            end
        end
        foreach(wait, tasks)
    end
    P
end
