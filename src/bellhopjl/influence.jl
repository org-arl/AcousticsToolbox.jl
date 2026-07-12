# SPDX-License-Identifier: GPL-3.0-or-later

# Geometric beam influence in Cartesian coordinates. Ports InfluenceGeoHatCart,
# InfluenceGeoGaussianCart, ApplyContribution, IncPhaseIfCaustic, IsAtCaustic
# and FinalPhase (influence.f90), and AddArr (ArrMod.f90).

"AD-safe equivalent of Fortran SPACING for the duplicate-point tests."
@inline _spacing(x) = 2.220446049250313e-16 * max(abs(x), one(x))

# IsAtCaustic with qleq0 = .TRUE. (the variant used by all geometric beams)
@inline _is_at_caustic(q, qold) =
    (q <= 0 && qold > 0) || (q >= 0 && qold < 0)

########################### sinks ###########################

abstract type InfluenceSink end

"Accumulates complex pressure (coherent) or intensity (in-/semicoherent)."
struct PressureSink{T<:Real} <: InfluenceSink
    U::Matrix{Complex{T}}     # (nz, nr)
    mode::RunMode
end

"One discrete arrival. Amplitude/phase/delay as stored by `AddArr` (ArrMod.f90)."
struct Arrival{T<:Real}
    amp::T
    phase::T
    delay::Complex{T}
    src_angle::T      # take-off angle [rad, positive down]
    rcv_angle::T      # receiver angle [rad, positive down]
    ntop::Int
    nbot::Int
end

"Accumulates discrete arrivals per receiver."
struct ArrivalSink{T<:Real} <: InfluenceSink
    arr::Matrix{Vector{Arrival{T}}}     # (nz, nr)
end

"""
Add an arrival, merging with the previous arrival at this receiver when both
delay and phase agree within `PhaseTol = 0.05` (ports `AddArr`, ArrMod.f90;
the amplitude-weighted merge keeps the previous phase, as the Fortran does).
"""
function add_arr!(v::Vector{Arrival{T}}, ω, amp, phase, delay, srcang, rcvang,
                  ntop, nbot) where {T}
    PhaseTol = T(0.05)
    if !isempty(v)
        b = v[end]
        if ω * abs(delay - b.delay) < PhaseTol && abs(b.phase - phase) < PhaseTol
            amptot = b.amp + amp
            w1 = b.amp / amptot
            w2 = amp / amptot
            v[end] = Arrival{T}(amptot, b.phase, w1 * b.delay + w2 * delay,
                                w1 * b.src_angle + w2 * srcang,
                                w1 * b.rcv_angle + w2 * rcvang, b.ntop, b.nbot)
            return v
        end
    end
    push!(v, Arrival{T}(amp, phase, delay, srcang, rcvang, ntop, nbot))
    v
end

# ApplyContribution (influence.f90)
@inline function apply_contribution!(sink::PressureSink{T}, iz, ir, amp, w,
                                     cnst, phaseInt, delay, ω, gaussian,
                                     srcang, rcvang, ntop, nbot) where {T}
    if sink.mode == COHERENT
        sink.U[iz, ir] += amp * exp(-im * (ω * delay - phaseInt))
    else  # incoherent / semicoherent: accumulate intensity
        s = (cnst * exp(imag(ω * delay)))^2 * w
        gaussian && (s *= sqrt(2 * T(π)))
        sink.U[iz, ir] += s
    end
    nothing
end

@inline function apply_contribution!(sink::ArrivalSink{T}, iz, ir, amp, w,
                                     cnst, phaseInt, delay, ω, gaussian,
                                     srcang, rcvang, ntop, nbot) where {T}
    add_arr!(sink.arr[iz, ir], ω, amp, phaseInt, delay, srcang, rcvang, ntop, nbot)
    nothing
end

########################### influence ###########################

"""
    influence_geo_cart!(sink, ray, α, Δα, ω, freq, rxr, rxz; gaussian=false)

Contribution of one beam to all receivers on the rectilinear grid
`rxr × rxz` (ranges ascending). Ports `InfluenceGeoHatCart` (`gaussian=false`)
and `InfluenceGeoGaussianCart` (`gaussian=true`) from influence.f90, including
the monotone receiver-pointer walk, the KMAH caustic phase tracking, and the
A-New-BellHope shallow-angle threshold change (0.50001) for Gaussian beams.
"""
function influence_geo_cart!(sink::InfluenceSink, ray::Vector{RayPt{T}}, α, Δα,
                             ω, freq, rxr::AbstractVector,
                             rxz::AbstractVector;
                             gaussian::Bool=false) where {T}
    nsteps = length(ray)
    nsteps < 2 && return nothing
    nrr = length(rxr)
    hugeT = T(Inf)

    q0 = ray[1].c / Δα                  # reference for J = q0 / q
    phase = zero(T)
    qOld = ray[1].q[1]
    rA = ray[1].x[1]

    # index of first receiver to the right of rA
    ir = searchsortedlast(rxr, rA) + 1
    ir > nrr && (ir = nrr)              # (Fortran MINLOC edge case; guarded)
    ray[1].t[1] < 0 && ir > 1 && (ir -= 1)   # left-travelling ray

    # sqrt(2π) represents a sum of Gaussians in free space; point source
    ratio1 = sqrt(abs(cos(α)))
    gaussian && (ratio1 /= sqrt(2 * T(π)))

    BeamWindow = 4      # Gaussian beam window: kills beams outside e^(-0.5 w²)
    srcang = T(α)

    @inbounds for is in 2:nsteps
        rB = ray[is].x[1]
        x_ray = ray[is-1].x

        rayt = ray[is].x - x_ray
        rlen = norm(rayt)
        rlen < 1.0e3 * _spacing(ray[is].x[1]) && continue   # duplicate point
        rayt = rayt / rlen
        rayn = SVector(-rayt[2], rayt[1])
        rcvang = atan(rayt[2], rayt[1])

        dqds = ray[is].q[1] - ray[is-1].q[1]
        dtauds = ray[is].τ - ray[is-1].τ

        qseg = ray[is-1].q[1]
        _is_at_caustic(qseg, qOld) && (phase += T(π) / 2)
        qOld = qseg

        λ = ray[is-1].c / freq
        if gaussian
            sigma = max(abs(ray[is-1].q[1]), abs(ray[is].q[1])) / q0 / abs(rayt[1])
            sigma = max(sigma, min(T(0.2) * freq * real(ray[is].τ), T(π) * λ))
            radiusMax = BeamWindow * sigma
            shallow = abs(rayt[1]) > T(0.50001)
        else
            radiusMax = max(abs(ray[is-1].q[1]), abs(ray[is].q[1])) / q0 / abs(rayt[1])
            shallow = abs(rayt[1]) > T(0.5)
        end

        # depth limits of beam
        if shallow
            zmin = min(ray[is-1].x[2], ray[is].x[2]) - radiusMax
            zmax = max(ray[is-1].x[2], ray[is].x[2]) + radiusMax
        else
            zmin = -hugeT
            zmax = +hugeT
        end

        while true
            # is Rr(ir) contained in [rA, rB)? then compute beam influence
            if min(rA, rB) <= rxr[ir] < max(rA, rB)
                for iz in eachindex(rxz)
                    x_rcvr = SVector(T(rxr[ir]), T(rxz[iz]))
                    (zmin <= x_rcvr[2] <= zmax) || continue

                    s = dot(x_rcvr - x_ray, rayt) / rlen  # proportional distance along ray
                    n = abs(dot(x_rcvr - x_ray, rayn))    # normal distance to ray
                    q = ray[is-1].q[1] + s * dqds         # interpolated amplitude

                    if gaussian
                        sigma = abs(q / q0)
                        sigma = max(sigma,
                                    min(T(0.2) * freq * real(ray[is].τ), T(π) * λ))
                        if n < BeamWindow * sigma
                            A = abs(q0 / q)
                            delay = ray[is-1].τ + s * dtauds
                            cnst = ratio1 * sqrt(ray[is].c / abs(q)) * ray[is].amp
                            w = exp(-(n / sigma)^2 / 2) / (sigma * A)
                            amp = cnst * w
                            # OALIB 2024 InfluenceGeoGaussianCart (influence.f90):
                            # ray phase read from the previous point (consistent
                            # with hat) and the caustic π/2 is additive (the 2022
                            # code discarded the accumulated phase)
                            phaseInt = ray[is-1].phase + phase
                            _is_at_caustic(q, qOld) && (phaseInt += T(π) / 2)
                            apply_contribution!(sink, iz, ir, amp, w, cnst,
                                                phaseInt, delay, ω, true,
                                                srcang, rcvang,
                                                ray[is].ntop, ray[is].nbot)
                        end
                    else
                        radius = abs(q / q0)              # beam radius
                        if n < radius
                            delay = ray[is-1].τ + s * dtauds
                            cnst = ratio1 * sqrt(ray[is].c / abs(q)) * ray[is].amp
                            w = (radius - n) / radius     # hat function
                            amp = cnst * w
                            # OALIB 2024 InfluenceGeoHatCart (influence.f90):
                            # caustic π/2 is additive, no longer discarding the
                            # accumulated ray phase
                            phaseInt = ray[is-1].phase + phase
                            _is_at_caustic(q, qOld) && (phaseInt += T(π) / 2)
                            apply_contribution!(sink, iz, ir, amp, w, cnst,
                                                phaseInt, delay, ω, false,
                                                srcang, rcvang,
                                                ray[is].ntop, ray[is].nbot)
                        end
                    end
                end
            end

            # bump receiver index towards rB
            if rxr[ir] < rB
                ir >= nrr && break
                rxr[ir + 1] >= rB && break
                ir += 1
            else
                ir <= 1 && break
                rxr[ir - 1] <= rB && break
                ir -= 1
            end
        end
        rA = rB
    end
    nothing
end
