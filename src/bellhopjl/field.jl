# SPDX-License-Identifier: GPL-3.0-or-later

# Beam loop, pressure assembly and arrivals. Ports the beam loop of
# BellhopCore and ScalePressure (influence.f90).

"Uniform fan of take-off angles (BellhopCore / ReadRayElevationAngles)."
function takeoff_angles(::Type{T}, beam::BeamParams) where {T}
    n = beam.nbeams
    n == 1 && return [T((beam.αmin + beam.αmax) / 2)], one(T)
    Δα = (T(beam.αmax) - T(beam.αmin)) / (n - 1)
    [T(beam.αmin) + (i - 1) * Δα for i in 1:n], Δα
end

"""
Automatic beam count for TL runs. Ports the Nalpha auto-selection of
`ReadRayElevationAngles` (angleMod.f90): phase-based limit with c0 = 1500 and
a beam-thinness limit w.r.t. water depth.
"""
function auto_nbeams(freq, rmax, depth)
    n = max(floor(Int, 0.3 * rmax * freq / 1500), 300)
    max(floor(Int, π / atan(depth / (10 * rmax))), n)
end

"Semi-coherent source amplitude (Lloyd mirror), BellhopCore ('S' runs)."
semicoherent_amp0(ω, c, zs, α) = sqrt(2) * abs(sin(ω / c * zs * sin(α)))

"Partition `1:n` into at most `k` contiguous, near-equal chunks."
_chunks(n, k) = [(1 + div((j - 1) * n, k)):div(j * n, k) for j in 1:min(k, n)]

"Preallocated ray-history buffer, reused across a beam loop."
function _ray_buffer(::Type{T}) where {T}
    hist = Vector{RayPt{T}}()
    sizehint!(hist, 2000)
    hist
end

"""
    pressure_field(env, zs, rxr, rxz, beam, mode; ntasks=1) -> Matrix{Complex} (nz × nr)

Trace all beams and assemble the receiver-grid pressure, including the final
`ScalePressure` step. `rxr` ascending ranges [m], `rxz` depths [m, +down].
With `ntasks > 1` the beam fan is split into contiguous chunks traced on
separate threads, each with its own accumulator; the per-chunk fields are
summed in chunk order, so results are reproducible for a given `ntasks`.
"""
function pressure_field(env::Env2D, zs, rxr, rxz, beam::BeamParams{TB},
                        mode::RunMode; ntasks::Int=1) where {TB}
    T = promote_type(env_eltype(env), typeof(float(zs)), TB)
    ω = 2π * env.freq
    αs, Δα = takeoff_angles(T, beam)
    e0, _ = evaluate(env.ssp, T(zs), one(T), 1)
    trace1! = (sink, hist, α) -> begin
        amp0 = mode == SEMICOHERENT ? semicoherent_amp0(ω, e0.c, T(zs), α) : one(T)
        ray = trace_ray!(hist, env, α, zs, beam; amp0)
        influence_geo_cart!(sink, ray, α, Δα, ω, env.freq, rxr, rxz;
                            gaussian=beam.gaussian)
    end
    U = zeros(Complex{T}, length(rxz), length(rxr))
    if ntasks <= 1 || length(αs) < 2 * ntasks
        sink = PressureSink{T}(U, mode)
        hist = _ray_buffer(T)
        for α in αs
            trace1!(sink, hist, α)
        end
    else
        tasks = map(_chunks(length(αs), ntasks)) do idxs
            Threads.@spawn begin
                local csink = PressureSink{T}(zeros(Complex{T}, length(rxz), length(rxr)), mode)
                local chist = _ray_buffer(T)
                for i in idxs
                    trace1!(csink, chist, αs[i])
                end
                csink.U
            end
        end
        for t in tasks              # fixed order ⇒ deterministic reduction
            U .+= fetch(t)
        end
    end
    scale_pressure!(U, mode, rxr)
end

"""
Final scaling of the pressure field. Ports `ScalePressure` (influence.f90) for
geometric beams (const = −1) and a point source: incoherent runs convert
intensity to pressure, then each column is scaled by −1/√r (0 at r = 0).
"""
function scale_pressure!(U::Matrix{Complex{T}}, mode::RunMode, rxr) where {T}
    if mode != COHERENT
        @. U = sqrt(complex(abs(real(U))))    # intensity → pressure
    end
    for (ir, r) in pairs(rxr)
        factor = iszero(r) ? zero(T) : -one(T) / sqrt(abs(T(r)))
        @views U[:, ir] .*= factor
    end
    U
end

"""
    compute_arrivals(env, zs, rxr, rxz, beam; ntasks=1) -> Matrix{Vector{Arrival}}

Arrivals at each receiver. Cylindrical-spreading factor 1/√r is applied to
amplitudes, as in `WriteArrivalsASCII` (ArrMod.f90). With `ntasks > 1` the
beam fan is split into contiguous chunks traced on separate threads, each
with its own sink; the per-chunk arrival lists are then replayed through
`add_arr!` in chunk (= beam) order, preserving the serial merge semantics.
"""
function compute_arrivals(env::Env2D, zs, rxr, rxz,
                          beam::BeamParams{TB}; ntasks::Int=1) where {TB}
    T = promote_type(env_eltype(env), typeof(float(zs)), TB)
    ω = 2π * env.freq
    αs, Δα = takeoff_angles(T, beam)
    arr = [Arrival{T}[] for _ in eachindex(rxz), _ in eachindex(rxr)]
    if ntasks <= 1 || length(αs) < 2 * ntasks
        sink = ArrivalSink{T}(arr)
        hist = _ray_buffer(T)
        for α in αs
            ray = trace_ray!(hist, env, α, zs, beam)
            influence_geo_cart!(sink, ray, α, Δα, ω, env.freq, rxr, rxz;
                                gaussian=beam.gaussian)
        end
    else
        tasks = map(_chunks(length(αs), ntasks)) do idxs
            Threads.@spawn begin
                local csink = ArrivalSink{T}([Arrival{T}[] for _ in eachindex(rxz),
                                                               _ in eachindex(rxr)])
                local chist = _ray_buffer(T)
                for i in idxs
                    local ray = trace_ray!(chist, env, αs[i], zs, beam)
                    influence_geo_cart!(csink, ray, αs[i], Δα, ω, env.freq,
                                        rxr, rxz; gaussian=beam.gaussian)
                end
                csink.arr
            end
        end
        for t in tasks              # replay in beam order (serial merge semantics)
            carr = fetch(t)
            for k in eachindex(arr)
                for a in carr[k]
                    add_arr!(arr[k], ω, a.amp, a.phase, a.delay, a.src_angle,
                             a.rcv_angle, a.ntop, a.nbot)
                end
            end
        end
    end
    for ir in eachindex(rxr), iz in eachindex(rxz)
        r = rxr[ir]
        factor = iszero(r) ? T(1e5) : 1 / sqrt(abs(T(r)))
        v = arr[iz, ir]
        for (k, a) in pairs(v)
            v[k] = Arrival{T}(factor * a.amp, a.phase, a.delay, a.src_angle,
                              a.rcv_angle, a.ntop, a.nbot)
        end
    end
    arr
end
