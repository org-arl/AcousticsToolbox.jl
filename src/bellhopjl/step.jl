# SPDX-License-Identifier: GPL-3.0-or-later

# Bellhop's modified midpoint stepper. Ports Step2D, ReduceStep2D and
# StepToBdry2D (Step.f90, A-New-BellHope 2022_4 layout including the
# boundary-snapping overhaul described in its README).

"""
    reduce_step2d(env, x0, urayt, iSegz0, topg, botg, box, deltas, h) -> (h, small)

Reduce trial step `h` along the unscaled tangent `urayt` so the step lands on
the nearest of: SSP depth-interface, ray box, top/bottom boundary, or
boundary-segment range interval. Ports `ReduceStep2D` (Step.f90).
"""
function reduce_step2d(env::Env2D, x0::SVector{2,T}, urayt::SVector{2,T},
                       iSegz0::Int, topg, botg, box, deltas, h) where {T}
    x = x0 + h * urayt
    hbig = T(Inf)

    # interface crossing in depth
    hInt = hbig
    zs = env.ssp.z
    if abs(urayt[2]) > eps(T)
        if zs[iSegz0] > x[2]
            hInt = (zs[iSegz0] - x0[2]) / urayt[2]
        elseif zs[iSegz0 + 1] < x[2]
            hInt = (zs[iSegz0 + 1] - x0[2]) / urayt[2]
        end
    end

    # ray box (centred at (0,0))
    hBoxr = hbig
    if abs(x[1]) > box.r
        hBoxr = (box.r - abs(x0[1])) / abs(urayt[1])
    end
    hBoxz = hbig
    if abs(x[2]) > box.z
        hBoxz = (box.z - abs(x0[2])) / abs(urayt[2])
    end

    # top crossing
    hTop = hbig
    d = x - topg.x
    if dot(topg.n, d) >= 0
        d0 = x0 - topg.x
        hTop = -dot(d0, topg.n) / dot(urayt, topg.n)
    end

    # bottom crossing
    hBot = hbig
    d = x - botg.x
    if dot(botg.n, d) >= 0
        d0 = x0 - botg.x
        hBot = -dot(d0, botg.n) / dot(urayt, botg.n)
    end

    # top or bottom segment crossing in range
    rSeg1 = max(topg.rseg[1], botg.rseg[1])
    rSeg2 = min(topg.rseg[2], botg.rseg[2])
    hSeg = hbig
    if abs(urayt[1]) > eps(T)
        if x[1] < rSeg1
            hSeg = -(x0[1] - rSeg1) / urayt[1]
        elseif x[1] > rSeg2
            hSeg = -(x0[1] - rSeg2) / urayt[1]
        end
    end

    h = min(h, hInt, hBoxr, hBoxz, hTop, hBot, hSeg)
    small = false
    if h < T(INFINITESIMAL_STEP_SIZE) * deltas
        h = T(INFINITESIMAL_STEP_SIZE) * deltas
        small = true
    end
    h, small
end

"""
    step_to_bdry2d(env, x0, urayt, iSegz0, topg, botg, box, deltas) ->
        (h, x2, topRefl, botRefl)

Take the blended tangent and find the minimum step size to put the ray exactly
on a boundary. Ports `StepToBdry2D` (Step.f90), including the exact snapping
to flat boundaries and the reflection flags.
"""
function step_to_bdry2d(env::Env2D, x0::SVector{2,T}, urayt::SVector{2,T},
                        iSegz0::Int, topg, botg, box, deltas) where {T}
    h = deltas
    x2 = x0 + h * urayt

    # interface crossing in depth
    zs = env.ssp.z
    if abs(urayt[2]) > eps(T)
        if zs[iSegz0] > x2[2]
            h = (zs[iSegz0] - x0[2]) / urayt[2]
            x2 = SVector(x0[1] + h * urayt[1], zs[iSegz0])
        elseif zs[iSegz0 + 1] < x2[2]
            h = (zs[iSegz0 + 1] - x0[2]) / urayt[2]
            x2 = SVector(x0[1] + h * urayt[1], zs[iSegz0 + 1])
        end
    end

    # ray box
    if abs(x2[1]) > box.r
        h = (box.r - abs(x0[1])) / abs(urayt[1])
        x2 = SVector(copysign(box.r, x0[1]), x0[2] + h * urayt[2])
    end
    if abs(x2[2]) > box.z
        h = (box.z - abs(x0[2])) / abs(urayt[2])
        x2 = SVector(x0[1] + h * urayt[1], copysign(box.z, x0[2]))
    end

    # top crossing
    topRefl = false
    d = x2 - topg.x
    if dot(topg.n, d) > -T(INFINITESIMAL_STEP_SIZE)
        d0 = x0 - topg.x
        h = -dot(d0, topg.n) / dot(urayt, topg.n)
        x2 = x0 + h * urayt
        if abs(topg.n[1]) < eps(T)                  # snap to exact depth if flat
            x2 = SVector(x2[1], topg.x[2])
        end
        topRefl = true
    end

    # bottom crossing
    botRefl = false
    d = x2 - botg.x
    if dot(botg.n, d) > -T(INFINITESIMAL_STEP_SIZE)
        d0 = x0 - botg.x
        h = -dot(d0, botg.n) / dot(urayt, botg.n)
        x2 = x0 + h * urayt
        if abs(botg.n[1]) < eps(T)
            x2 = SVector(x2[1], botg.x[2])
        end
        botRefl = true
        topRefl = false
    end

    # top or bottom segment crossing in range
    rSeg1 = max(topg.rseg[1], botg.rseg[1])
    rSeg2 = min(topg.rseg[2], botg.rseg[2])
    if abs(urayt[1]) > eps(T)
        if x2[1] < rSeg1
            h = -(x0[1] - rSeg1) / urayt[1]
            x2 = SVector(rSeg1, x0[2] + h * urayt[2])
            topRefl = false
            botRefl = false
        elseif x2[1] > rSeg2
            h = -(x0[1] - rSeg2) / urayt[1]
            x2 = SVector(rSeg2, x0[2] + h * urayt[2])
            topRefl = false
            botRefl = false
        end
    end

    if h < T(INFINITESIMAL_STEP_SIZE) * deltas       # infinitesimal step
        h = T(INFINITESIMAL_STEP_SIZE) * deltas
        x2 = x0 + h * urayt
        # recheck reflection conditions
        topRefl = dot(topg.n, x2 - topg.x) > eps(T)
        botRefl = dot(botg.n, x2 - botg.x) > eps(T)
        botRefl && (topRefl = false)
    end
    h, x2, topRefl, botRefl
end

@inline cnn_csq(e::SSPEval, t::SVector{2}) = e.czz * t[1]^2
# (crr·t₂² − 2·crz·t₁·t₂ + czz·t₁² with crr = crz = 0 for range-independent SSP)

"""
    step2d(env, ray0, iSegz, topg, botg, box, deltas) ->
        (ray2, iSegz2, topRefl, botRefl, small)

One Bellhop integration step. Ports `Step2D` (Step.f90): reduced half-Euler to
the midpoint, second `reduce_step2d` with midpoint derivatives, weighted
update via `step_to_bdry2d`, then the SSP interface jump condition on `p`.
"""
function step2d(env::Env2D, ray0::RayPt{T}, iSegz::Int, topg, botg, box,
                deltas) where {T}
    # phase 1 (Euler half step)
    e0, iSegz = evalssp(env, ray0.x, ray0.t, iSegz)
    c0 = e0.c
    csq0 = c0 * c0
    cnn0_csq0 = cnn_csq(e0, ray0.t)
    iSegz0 = iSegz
    gradc0 = SVector(zero(e0.cz), e0.cz)

    h = deltas
    urayt0 = c0 * ray0.t
    h, small1 = reduce_step2d(env, ray0.x, urayt0, iSegz0, topg, botg, box, deltas, h)
    halfh = h / 2

    x1 = ray0.x + halfh * urayt0
    t1 = ray0.t - halfh * gradc0 / csq0
    p1 = ray0.p - halfh * cnn0_csq0 * ray0.q
    q1 = ray0.q + halfh * c0 * ray0.p

    # phase 2
    e1, iSegz = evalssp(env, x1, t1, iSegz)
    c1 = e1.c
    csq1 = c1 * c1
    cnn1_csq1 = cnn_csq(e1, t1)
    gradc1 = SVector(zero(e1.cz), e1.cz)

    urayt1 = c1 * t1
    h, small2 = reduce_step2d(env, ray0.x, urayt1, iSegz0, topg, botg, box, deltas, h)

    # blend of f' based on proportion of a full step used
    w1 = h / (2 * halfh)
    w0 = 1 - w1
    urayt2   =  w0 * urayt0 + w1 * urayt1
    unitdt   = -w0 * gradc0 / csq0 - w1 * gradc1 / csq1
    unitdp   = -w0 * cnn0_csq0 * ray0.q - w1 * cnn1_csq1 * q1
    unitdq   =  w0 * c0 * ray0.p + w1 * c1 * p1
    unitdtau =  w0 / complex(c0, e0.cimag) + w1 / complex(c1, e1.cimag)

    h, x2, topRefl, botRefl =
        step_to_bdry2d(env, ray0.x, urayt2, iSegz0, topg, botg, box, deltas)
    t2 = ray0.t + h * unitdt
    p2 = ray0.p + h * unitdp
    q2 = ray0.q + h * unitdq
    τ2 = ray0.τ + h * unitdtau

    e2, iSegz = evalssp(env, x2, t2, iSegz)

    # if we crossed an interface, apply jump condition on p
    # (OALIB 2024 Step2D, Step.f90: only for SSPs with a discontinuous first
    # derivative — 'N'/'C'; spline/PCHIP/analytic profiles are C¹ and get none)
    if has_gradcjump(env.ssp) && iSegz != iSegz0 && !iszero(t2[2])
        gradcjump = SVector(zero(e2.cz), e2.cz - e0.cz)
        ray2n = SVector(-t2[2], t2[1])
        cnjump = dot(gradcjump, ray2n)
        csjump = dot(gradcjump, t2)
        RM = t2[1] / t2[2]                # crossing in depth
        RN = RM * (2 * cnjump - RM * csjump) / e2.c
        p2 = p2 - q2 * RN
    end

    ray2 = RayPt{T}(x2, t2, p2, q2, τ2, e2.c, ray0.amp, ray0.phase,
                    ray0.ntop, ray0.nbot)
    ray2, iSegz, topRefl, botRefl, small1 | small2
end
