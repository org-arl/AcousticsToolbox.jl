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

"""
    pressure_field(env, zs, rxr, rxz, beam, mode) -> Matrix{Complex} (nz × nr)

Trace all beams and assemble the receiver-grid pressure, including the final
`ScalePressure` step. `rxr` ascending ranges [m], `rxz` depths [m, +down].
"""
function pressure_field(env::Env2D, zs, rxr, rxz, beam::BeamParams{TB},
                        mode::RunMode) where {TB}
    T = promote_type(env_eltype(env), typeof(float(zs)), TB)
    ω = 2π * env.freq
    αs, Δα = takeoff_angles(T, beam)
    sink = PressureSink{T}(zeros(Complex{T}, length(rxz), length(rxr)), mode)
    e0, _ = evaluate(env.ssp, T(zs), one(T), 1)
    for α in αs
        amp0 = mode == SEMICOHERENT ? semicoherent_amp0(ω, e0.c, T(zs), α) : one(T)
        ray = trace_ray(env, α, zs, beam; amp0)
        influence_geo_cart!(sink, ray, α, Δα, ω, env.freq, rxr, rxz;
                            gaussian=beam.gaussian)
    end
    scale_pressure!(sink.U, mode, rxr)
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
    compute_arrivals(env, zs, rxr, rxz, beam) -> Matrix{Vector{Arrival}}

Arrivals at each receiver. Cylindrical-spreading factor 1/√r is applied to
amplitudes, as in `WriteArrivalsASCII` (ArrMod.f90).
"""
function compute_arrivals(env::Env2D, zs, rxr, rxz,
                          beam::BeamParams{TB}) where {TB}
    T = promote_type(env_eltype(env), typeof(float(zs)), TB)
    ω = 2π * env.freq
    αs, Δα = takeoff_angles(T, beam)
    arr = [Arrival{T}[] for _ in eachindex(rxz), _ in eachindex(rxr)]
    sink = ArrivalSink{T}(arr)
    for α in αs
        ray = trace_ray(env, α, zs, beam)
        influence_geo_cart!(sink, ray, α, Δα, ω, env.freq, rxr, rxz;
                            gaussian=beam.gaussian)
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
