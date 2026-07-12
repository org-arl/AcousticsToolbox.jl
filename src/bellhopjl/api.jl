# SPDX-License-Identifier: GPL-3.0-or-later

# UnderwaterAcoustics.jl integration. This is the only file that touches
# UnderwaterAcoustics.jl types; it mirrors the environment handling (and unit
# conversions) of the Fortran Bellhop wrapper (src/bellhop.jl and
# src/common.jl _write_env) so that BellhopJL is a drop-in alternative.

using UnderwaterAcoustics
using UnderwaterAcoustics: AbstractRayPropagationModel, AbstractAcousticSource,
                           AbstractAcousticReceiver, AcousticReceiverGrid2D,
                           RayArrival, SampledFieldX, SampledFieldZ,
                           RigidBoundary, PressureReleaseBoundary, FluidBoundary,
                           is_range_dependent, is_constant, value, in_units,
                           db2amp, in_dBperλ, CubicSpline, Linear, XYZ

"""
    BellhopJL(env; kwargs...)

Native-Julia port of the 2D Bellhop Gaussian-beam/ray tracer — a drop-in
alternative to the Fortran-wrapping `Bellhop` model in this package,
differentiable with ForwardDiff.

Supported keyword arguments:
- `nbeams`: number of beams to use (default: 0, auto)
- `min_angle`: minimum beam angle (default: -80°)
- `max_angle`: maximum beam angle (default: 80°)
- `beam_type`: `:geometric` (default) or `:gaussian`
"""
struct BellhopJL{T} <: AbstractRayPropagationModel
    env::T
    nbeams::Int
    min_angle::Float64      # [rad], UA convention (positive up)
    max_angle::Float64
    gaussian::Bool
    function BellhopJL(env, nbeams, min_angle, max_angle, beam_type)
        env.seabed isa FluidBoundary ||
            error("seabed must be a FluidBoundary")
        env.surface isa FluidBoundary ||
            error("surface must be a FluidBoundary")
        min_angle = in_units(u"rad", min_angle)
        max_angle = in_units(u"rad", max_angle)
        -π/2 ≤ min_angle ≤ π/2 || error("min_angle should be between -π/2 and π/2")
        -π/2 ≤ max_angle ≤ π/2 || error("max_angle should be between -π/2 and π/2")
        min_angle < max_angle || error("max_angle should be more than min_angle")
        beam_type ∈ (:geometric, :gaussian) || error("Unknown beam_type type")
        new{typeof(env)}(env, max(nbeams, 0), Float64(min_angle),
                         Float64(max_angle), beam_type === :gaussian)
    end
end

BellhopJL(env; nbeams=0, min_angle=-80°, max_angle=80°, beam_type=:geometric) =
    BellhopJL(env, nbeams, min_angle, max_angle, beam_type)

Base.show(io::IO, pm::BellhopJL) = print(io, "BellhopJL(⋯)")

### adapter helpers

function _check2d(tx, rxs)
    location(tx).x == 0 || error("Transmitter must be at (0, 0, z)")
    location(tx).y == 0 || error("2D model requires transmitter in the x-z plane")
    all(location(rx).x >= 0 for rx in rxs) || error("Receivers must be in the +x halfspace")
    all(location(rx).y == 0 for rx in rxs) || error("2D model requires receivers in the x-z plane")
    nothing
end

# nominal half-wavelength sampling for range-dependent boundary functions
# (mirrors _recommend_len in AcousticsToolbox.jl/src/common.jl)
function _recommend_len(x, f)
    λ = 1500.0 / f
    clamp(round(Int, 2x / λ) + 1, 25, 1000)
end

# UA FluidBoundary → Bellhop halfspace BC. Unit conversions copied from the
# wrapper's env writer: density ratio = ρ/1000 (g/cm³, water = 1), attenuation
# in dB/λ via in_dBperλ(δ).
function _make_bc(b::FluidBoundary, f)
    if isinf(b.c) && isinf(b.ρ)
        RigidBC()
    elseif b.c == 0 && b.ρ == 0
        VacuumBC()
    else
        HalfspaceBC(b.ρ / 1000, b.c, in_dBperλ(b.δ), f)
    end
end

# UA soundspeed field → (depth nodes, complex speed nodes, spline?)
# mirrors the SSP block of _write_env (common.jl)
function _make_ssp_nodes(ssp, waterdepth, f, fg)
    if is_constant(ssp)
        c = value(ssp)
        zn = [zero(waterdepth), waterdepth]
        cn = [crci(c, zero(c), f; volume=fg) for _ in 1:2]
        return zn, cn, false
    elseif ssp isa SampledFieldZ
        zrange = sort!(vcat(collect(ssp.zrange), -waterdepth); rev=true)
        zn = eltype(zrange)[]
        cn = []
        for z in zrange
            push!(zn, -z)
            push!(cn, crci(ssp(z), zero(ssp(z)), f; volume=fg))
            z == -waterdepth && break
        end
        return zn, [c for c in cn], ssp.interp === CubicSpline()
    else
        zn = collect(range(zero(waterdepth), waterdepth; length=_recommend_len(waterdepth, f)))
        cn = [crci(ssp(-z), zero(ssp(-z)), f; volume=fg) for z in zn]
        if !any(z -> z == waterdepth, zn)
            push!(zn, waterdepth)
            push!(cn, crci(ssp(-waterdepth), zero(ssp(-waterdepth)), f; volume=fg))
        end
        return zn, cn, false
    end
end

# UA bathymetry/altimetry → Boundary2D (mirrors _create_alt_bathy_file sampling)
function _make_boundary(data, func, maxr, f; top::Bool)
    if !is_range_dependent(data)
        return flat_boundary(func(data, (0.0, 0.0)); top)
    end
    x = data isa SampledFieldX ? collect(data.xrange) :
        collect(range(0.0, maxr; length=_recommend_len(maxr, f)))
    z = [func(data, (x1, 0.0)) for x1 in x]
    Boundary2D(x, z; top)
end

"""
Build the internal `Env2D` (+ source depth, beam parameters) from a UA
environment, source and receiver extent.
"""
function _make_env2d(pm::BellhopJL, tx, maxr; nbeams=0)
    env = pm.env
    f = frequency(tx)
    zs = -location(tx).z
    bathy = env.bathymetry
    waterdepth = maximum(bathy)
    fg = FrancoisGarrison(float(env.temperature), float(env.salinity),
                          float(env.pH), float(waterdepth / 2))

    zn, cn, spline = _make_ssp_nodes(env.soundspeed, waterdepth, f, fg)
    ssp = spline ? SplineSSP(zn, cn) : CLinearSSP(zn, cn)

    top = _make_boundary(env.altimetry, (q, p) -> -value(q, p), maxr, f; top=true)
    bot = _make_boundary(bathy, value, maxr, f; top=false)

    topbc = _make_bc(env.surface, f)
    botbc = _make_bc(env.seabed, f)

    env2d = Env2D(float(f), ssp, top, bot, topbc, botbc)

    # element type for beam parameters (must carry duals of bathymetry/range/freq)
    T = promote_type(eltype(ssp), typeof(float(waterdepth)), typeof(float(maxr)),
                     typeof(float(f)))
    n = nbeams > 0 ? nbeams :
        pm.nbeams > 0 ? pm.nbeams : auto_nbeams(f, maxr, waterdepth)
    beam = BeamParams(deltas=T(waterdepth / 10), box_r=T(1.01 * maxr),
                      box_z=T(1.01 * waterdepth), nbeams=n,
                      αmin=T(-pm.max_angle), αmax=T(-pm.min_angle),
                      gaussian=pm.gaussian)
    env2d, zs, beam
end

_runmode(mode) = mode === :coherent ? COHERENT :
                 mode === :incoherent ? INCOHERENT :
                 mode === :semicoherent ? SEMICOHERENT :
                 error("Unknown mode :" * string(mode))

### interface functions

function UnderwaterAcoustics.acoustic_field(pm::BellhopJL,
        tx::AbstractAcousticSource, rxs::AcousticReceiverGrid2D; mode=:coherent)
    _check2d(tx, rxs)
    runmode = _runmode(mode)
    xrev = first(rxs.xrange) > last(rxs.xrange)
    rxr = xrev ? reverse(collect(rxs.xrange)) : collect(rxs.xrange)
    rxz = [-z for z in rxs.zrange]
    maxr = maximum(rxs.xrange)
    env2d, zs, beam = _make_env2d(pm, tx, maxr)
    U = pressure_field(env2d, zs, rxr, rxz, beam, runmode)   # (nz, nr)
    fld = -permutedims(U) .* db2amp(spl(tx))                 # (nr, nz)
    xrev ? reverse(fld; dims=1) : fld
end

function UnderwaterAcoustics.acoustic_field(pm::BellhopJL,
        tx::AbstractAcousticSource, rx::AbstractAcousticReceiver; mode=:coherent)
    _check2d(tx, (rx,))
    runmode = _runmode(mode)
    p = location(rx)
    env2d, zs, beam = _make_env2d(pm, tx, p.x)
    U = pressure_field(env2d, zs, [p.x], [-p.z], beam, runmode)
    -U[1, 1] * db2amp(spl(tx))
end

function UnderwaterAcoustics.arrivals(pm::BellhopJL,
        tx::AbstractAcousticSource, rx::AbstractAcousticReceiver; paths=true)
    _check2d(tx, (rx,))
    p = location(rx)
    nbeams = pm.nbeams
    nbeams == 0 &&
        (nbeams = round(Int, (pm.max_angle - pm.min_angle) / deg2rad(0.05)) + 1)
    env2d, zs, beam = _make_env2d(pm, tx, p.x; nbeams)
    f = env2d.freq
    arr = compute_arrivals(env2d, zs, [p.x], [-p.z], beam)[1, 1]
    out = map(arr) do a
        ϕ = a.amp * cis(a.phase) * exp(2π * f * imag(a.delay))
        path = if paths
            ray = trace_ray(env2d, a.src_angle, zs, beam)
            [(x=pt.x[1], y=zero(pt.x[1]), z=-pt.x[2]) for pt in ray]
        else
            missing
        end
        RayArrival(real(a.delay), complex(ϕ), a.ntop, a.nbot,
                   -a.src_angle, a.rcv_angle, path)
    end
    sort(out; by=a -> a.t)
end
