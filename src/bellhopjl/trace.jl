# SPDX-License-Identifier: GPL-3.0-or-later

# Ray tracing and boundary reflection. Ports TraceRay2D and Reflect2D
# (bellhop.f90).

"""
    reflect2d(env, ray, bc, tBdry, nBdry, κ, ω, iSegz, is_top) -> (RayPt, iSegz)

Specular reflection at a boundary with unit tangent `tBdry`, outward unit
normal `nBdry` and curvature `κ`. Ports `Reflect2D` (bellhop.f90) with the
default beam type `Beam%Type(3:3) = 'S'` (standard curvature condition — no
doubling/zeroing) and no beam shift.
"""
function reflect2d(env::Env2D, ray::RayPt{T}, bc::BoundaryCondition,
                   tBdry::SVector{2}, nBdry::SVector{2}, κ, ω, iSegz::Int,
                   is_top::Bool) where {T}
    Tg = dot(ray.t, tBdry)      # component of ray tangent along boundary
    Th = dot(ray.t, nBdry)      # component of ray tangent normal to boundary

    t′ = ray.t - 2 * Th * nBdry  # reflected (scaled) tangent

    # curvature correction, based on Muller, Geoph. J. R.A.S., 79 (1984)
    e, iSegz = evalssp(env, ray.x, ray.t, iSegz)   # just to get c
    c = e.c
    gradc = SVector(zero(e.cz), e.cz)

    rayt = c * ray.t                           # incident unit tangent
    rayn = SVector(-rayt[2], rayt[1])          # incident unit normal
    rayt_tilde = c * t′                        # reflected unit tangent
    rayn_tilde = -SVector(-rayt_tilde[2], rayt_tilde[1])  # reflected unit normal

    RN = 2 * κ / c^2 / Th                      # boundary curvature correction

    cnjump = -dot(gradc, rayn_tilde - rayn)
    csjump = -dot(gradc, rayt_tilde - rayt)

    if is_top
        cnjump = -cnjump   # (t,n) system of the top boundary has opposite sense
        RN = -RN
    end

    RM = Tg / Th           # tan(angle of incidence)
    RN = RN + RM * (2 * cnjump - RM * csjump) / c^2
    # Beam%Type(3:3) = 'S' (default): standard curvature condition (no scaling)

    p′ = ray.p + ray.q * RN
    q′ = ray.q

    aR, φR = reflection_coefficient(bc, Tg, Th, ω)
    amp′ = iszero(aR) ? zero(ray.amp) : aR * ray.amp
    phase′ = ray.phase + φR

    RayPt{T}(ray.x, t′, p′, q′, ray.τ, c, amp′, phase′,
             ray.ntop + (is_top ? 1 : 0), ray.nbot + (is_top ? 0 : 1)), iSegz
end

"""
    trace_ray(env, α, zs, beam; amp0=1) -> Vector{RayPt}

Trace the beam with take-off angle `α` [rad, positive down] from the source at
(0, zs). Ports `TraceRay2D` (bellhop.f90). At each boundary reflection two
coincident points are stored (incoming and outgoing tangents), as in the
Fortran.
"""
function trace_ray(env::Env2D, α, zs, beam::BeamParams; amp0=1)
    T = promote_type(env_eltype(env), typeof(float(α)), typeof(float(zs)),
                     typeof(beam.deltas), typeof(float(amp0)))
    hist = Vector{RayPt{T}}()
    sizehint!(hist, 2000)
    trace_ray!(hist, env, α, zs, beam; amp0)
end

"""
    trace_ray!(hist, env, α, zs, beam; amp0=1) -> hist

In-place variant of [`trace_ray`](@ref): empties `hist` and fills it with the
traced ray, so one buffer can be reused across the beam loop.
"""
function trace_ray!(hist::Vector{RayPt{T}}, env::Env2D, α, zs,
                    beam::BeamParams; amp0=1) where {T}
    empty!(hist)
    ω = 2π * env.freq
    deltas = T(beam.deltas)
    box = (r = T(beam.box_r), z = T(beam.box_z))

    iSegz = 1
    xs = SVector{2,T}(zero(T), T(zs))
    tinit = SVector{2,T}(cos(α), sin(α))
    e0, iSegz = evalssp(env, xs, tinit, iSegz)

    ray = RayPt{T}(xs, tinit / e0.c, SVector{2,T}(1, 0), SVector{2,T}(0, 0),
                   zero(Complex{T}), e0.c, T(amp0), zero(T), 0, 0)
    # q initialised to (0,0): TraceRay2D sets q = 0 for geometric beam runs
    # ('G'); the (p₂,q₂) plane-wave pair is decoupled from (p₁,q₁) and unused
    # by the geometric influence functions, so this also matches 'B' runs.
    push!(hist, ray)

    iTop = get_seg(env.top, xs[1], ray.t[1])
    iBot = get_seg(env.bot, xs[1], ray.t[1])
    topg = segment_geom(env.top, iTop)
    botg = segment_geom(env.bot, iBot)

    distBegTop, distBegBot = distances(ray.x, topg, botg)
    if distBegTop <= 0 || distBegBot <= 0
        return hist   # source on or outside the boundaries
    end

    while true
        ray, iSegz, topRefl, botRefl, _ =
            step2d(env, ray, iSegz, topg, botg, box, deltas)
        push!(hist, ray)

        # new altimetry segment?
        r = ray.x[1]
        if r < topg.rseg[1] || (r == topg.rseg[1] && ray.t[1] < 0) ||
           r > topg.rseg[2] || (r == topg.rseg[2] && ray.t[1] >= 0)
            iTop = get_seg(env.top, r, ray.t[1])
            topg = segment_geom(env.top, iTop)
        end
        # new bathymetry segment?
        if r < botg.rseg[1] || (r == botg.rseg[1] && ray.t[1] < 0) ||
           r > botg.rseg[2] || (r == botg.rseg[2] && ray.t[1] >= 0)
            iBot = get_seg(env.bot, r, ray.t[1])
            botg = segment_geom(env.bot, iBot)
        end

        distEndTop, distEndBot = distances(ray.x, topg, botg)

        if topRefl
            ray, iSegz = reflect2d(env, ray, env.topbc, topg.t, topg.n, topg.κ,
                                   ω, iSegz, true)
            push!(hist, ray)
            distEndTop, distEndBot = distances(ray.x, topg, botg)
        elseif botRefl
            ray, iSegz = reflect2d(env, ray, env.botbc, botg.t, botg.n, botg.κ,
                                   ω, iSegz, false)
            push!(hist, ray)
            distEndTop, distEndBot = distances(ray.x, topg, botg)
        end

        # has the ray left the box, lost its energy, or escaped the boundaries?
        if abs(ray.x[1]) >= box.r || abs(ray.x[2]) >= box.z ||
           ray.amp < T(0.005) ||
           (distBegTop < 0 && distEndTop < 0) ||
           (distBegBot < 0 && distEndBot < 0) ||
           length(hist) >= MAX_N - 3
            break
        end

        distBegTop = distEndTop
        distBegBot = distEndBot
    end
    hist
end
