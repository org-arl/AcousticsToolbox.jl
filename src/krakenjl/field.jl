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
    iz = [clamp(searchsortedlast(z, zt), 1, n - 1) for zt in ztab]
    w = [(ztab[k] - z[iz[k]]) / (z[iz[k]+1] - z[iz[k]]) for k in eachindex(iz)]
    iz, w
end

"Interpolate mode shapes onto depths `ztab` (Vector, kraken(c).f90 PhiTab)."
function tabulate_modes(res::ModeResult{T}, ztab::AbstractVector) where {T}
    iz, wts = weight(res.z, ztab)
    M = length(res.k)
    TT = promote_type(T, typeof(float(one(eltype(ztab)))))
    phitab = Matrix{Complex{TT}}(undef, length(ztab), M)
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
function evaluate_field(k::Vector{<:Complex}, phiS::Vector{<:Complex},
                        phiR::Matrix{<:Complex}, rxr::AbstractVector;
                        incoherent::Bool=false, threads::Int=1)
    T = promote_type(real(eltype(k)), real(eltype(phiS)), real(eltype(phiR)),
                     typeof(float(one(eltype(rxr)))))
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
            P[iz, ir] = abs(_value(r)) > floatmin(_valuetype(T)) ? s / sqrt(r) : s
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
