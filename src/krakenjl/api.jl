# SPDX-License-Identifier: GPL-3.0-or-later

# UnderwaterAcoustics.jl integration. This is the only file that touches
# UnderwaterAcoustics.jl types; it mirrors the environment handling (and unit
# conversions) of the Fortran Kraken wrapper (src/kraken.jl and
# src/common.jl _write_env) so that KrakenJL is a drop-in alternative.

using UnderwaterAcoustics
using UnderwaterAcoustics: AbstractModePropagationModel, AbstractAcousticSource,
                           AbstractAcousticReceiver, AcousticReceiverGrid2D,
                           ModeArrival, SampledFieldZ, SampledFieldXZ,
                           RigidBoundary, PressureReleaseBoundary, FluidBoundary,
                           ElasticBoundary, MultilayerElasticBoundary,
                           is_range_dependent, is_constant, value, in_units,
                           db2amp, in_dBperλ, CubicSpline

"""
    KrakenJL(env; kwargs...)

Native-Julia port of the KRAKEN / KRAKENC normal-mode models — a drop-in
alternative to the Fortran-wrapping `Kraken` model in this package.

Supported keyword arguments:
- `nmodes`: number of modes to use (default: 9999)
- `mesh_density`: number of mesh points per wavelength (default: 0, 0=auto)
- `clow`: lower limit of phase speed (default: 1300, 0=auto)
- `chigh`: upper limit of phase speed (default: 2500)
- `rmax`: largest range (in m) of interest (default: Inf, meaning 1% more
  than the furthest receiver)
- `complex_solver`: use the KRAKENC complex solver (default: true)
- `robust`: use robust (but slow) root finder restarts (default: false)
- `threads`: number of threads for mode/field computation
  (default: `Threads.nthreads()`; 1 runs the serial code path)
"""
struct KrakenJL{T} <: AbstractModePropagationModel
    env::T
    nmodes::Int
    mesh_density::Float64
    clow::Float64
    chigh::Float64
    rmax::Float64
    complex_solver::Bool
    robust::Bool
    threads::Int
    function KrakenJL(env, nmodes, mesh_density, clow, chigh, rmax,
                      complex_solver, robust, threads)
        _check_krakenjl_env(env)
        nmodes ≥ 1 || error("number of modes should be positive")
        mesh_density ≥ 0 || error("mesh density should be non-negative")
        clow ≥ 0 || error("clow should be non-negative")
        chigh > clow || error("chigh should be more than clow")
        rmax ≥ 0 || error("rmax should be non-negative")
        threads ≥ 1 || error("threads should be at least 1")
        new{typeof(env)}(env, nmodes, Float64(mesh_density), Float64(clow),
                         Float64(chigh), Float64(rmax), complex_solver, robust,
                         threads)
    end
end

KrakenJL(env; nmodes=9999, mesh_density=0, clow=1300.0, chigh=2500.0, rmax=Inf,
         complex_solver=true, robust=false, threads=Threads.nthreads()) =
    KrakenJL(env, nmodes, mesh_density, clow, chigh, rmax, complex_solver,
             robust, threads)

Base.show(io::IO, pm::KrakenJL) = print(io, "KrakenJL(⋯)")

### adapter helpers

# mirrors _check_env(::Type{Kraken}, env) in src/kraken.jl
function _check_krakenjl_env(env)
    env.seabed isa FluidBoundary || env.seabed isa ElasticBoundary ||
        env.seabed isa MultilayerElasticBoundary ||
        error("seabed must be a FluidBoundary, ElasticBoundary or MultilayerElasticBoundary")
    env.surface isa FluidBoundary || error("surface must be a FluidBoundary")
    is_range_dependent(env.soundspeed) && error("Range-dependent soundspeed not supported")
    is_range_dependent(env.altimetry) && error("Range-dependent altimetry not supported")
    is_range_dependent(env.bathymetry) && error("Range-dependent bathymetry not supported")
    nothing
end

# promoted scalar type for the solve, so ForwardDiff duals in the frequency,
# soundspeed, bathymetry or boundary parameters flow through
function _kj_eltype(env, f)
    T = promote_type(typeof(float(f)), typeof(float(maximum(env.bathymetry))),
                     typeof(float(maximum(env.soundspeed))), Float64)
    for b in (env.surface, env.seabed)
        if b isa FluidBoundary
            T = promote_type(T, typeof(float(b.c)), typeof(float(b.ρ)),
                             typeof(float(b.δ)))
        elseif b isa ElasticBoundary
            T = promote_type(T, typeof(float(b.cₚ)), typeof(float(b.cₛ)),
                             typeof(float(b.ρ)), typeof(float(b.δₚ)),
                             typeof(float(b.δₛ)))
        elseif b isa MultilayerElasticBoundary
            for l in b.layers
                T = promote_type(T, typeof(float(first(l.cₚ))),
                                 typeof(float(first(l.cₛ))),
                                 typeof(float(first(l.ρ))), typeof(float(l.h)))
            end
        end
    end
    T
end

# nominal half-wavelength sampling (mirrors _recommend_len, src/common.jl)
function _kj_recommend_len(x, f)
    λ = 1500.0 / f
    clamp(round(Int, 2x / λ) + 1, 25, 1000)
end

# UA boundary → Halfspace. Unit conversions mirror _write_env (src/common.jl):
# density ratio ρ/1000 (g/cm³), attenuation in dB/λ, Francois-Garrison volume
# attenuation applied through crci as the Fortran CRCI does for AttenUnit 'WF'
function _kj_halfspace(::Type{T}, b, f, fg) where {T}
    if b === RigidBoundary
        Halfspace{T}(:rigid)
    elseif b === PressureReleaseBoundary
        Halfspace{T}(:vacuum)
    elseif b isa ElasticBoundary
        if b.cₛ == 0    # cₛ = 0 elastic reduces to fluid (mirrors _write_env)
            Halfspace{T}(:acoustoelastic,
                         crci(float(b.cₚ), in_dBperλ(b.δₚ), f; volume=fg),
                         zero(Complex{T}), b.ρ / 1000)
        else
            Halfspace{T}(:acoustoelastic,
                         crci(float(b.cₚ), in_dBperλ(b.δₚ), f; volume=fg),
                         crci(float(b.cₛ), in_dBperλ(b.δₛ), f; volume=fg),
                         b.ρ / 1000)
        end
    else                # FluidBoundary
        Halfspace{T}(:acoustoelastic,
                     crci(float(b.c), in_dBperλ(b.δ), f; volume=fg),
                     zero(Complex{T}), b.ρ / 1000)
    end
end

"""
Build the internal `Problem` (+ base mesh counts) from a UA environment,
source and receiver extent. Mirrors the Kraken blocks of `_write_env`
(src/common.jl) and the auto-mesh rule of ReadEnvironmentMod.f90.
"""
function _make_problem(pm::KrakenJL, tx, maxr)
    env = pm.env
    f = float(frequency(tx))
    T = _kj_eltype(env, f)
    ssp = env.soundspeed
    waterdepth = T(maximum(env.bathymetry))
    fg = FrancoisGarrison(float(env.temperature), float(env.salinity),
                          float(env.pH), float(waterdepth / 2))

    media = Medium{T}[]
    sigma = T[float(env.surface.σ)]

    # water column (SSP block of _write_env; water density ratio = 1, water
    # attenuation = volume attenuation only)
    spline = ssp isa SampledFieldZ && ssp.interp === CubicSpline()
    zn = T[]
    cn = Complex{T}[]
    if is_constant(ssp)
        c = float(value(ssp))
        zn = [zero(T), waterdepth]
        cn = Complex{T}[crci(c, 0.0, f; volume=fg) for _ in 1:2]
    elseif ssp isa SampledFieldZ
        zrange = sort!(vcat(collect(ssp.zrange), -waterdepth); rev=true)
        for z in zrange
            push!(zn, -z)
            push!(cn, crci(float(ssp(z)), 0.0, f; volume=fg))
            z == -waterdepth && break
        end
    else
        for d in range(zero(T), waterdepth;
                       length=_kj_recommend_len(_value(waterdepth), _value(f)))
            push!(zn, d)
            push!(cn, crci(float(ssp(-d)), 0.0, f; volume=fg))
        end
        if floor(waterdepth) != waterdepth
            push!(zn, waterdepth)
            push!(cn, crci(float(ssp(-waterdepth)), 0.0, f; volume=fg))
        end
    end
    λ = maximum(ssp) / f
    ng_water = round(Int, _value(waterdepth / λ) * pm.mesh_density)
    push!(media, Medium{T}(0.0, waterdepth, zn, cn, zeros(Complex{T}, length(zn)),
                           ones(T, length(zn)), spline, ng_water))

    # sediment layers of a multilayer elastic seabed (mirrors _write_env)
    depth = waterdepth
    seabed = env.seabed
    if seabed isa MultilayerElasticBoundary
        for l in seabed.layers[1:end-1]
            ρ₁, ρ₂ = float(first(l.ρ)), float(last(l.ρ))
            cₚ₁, cₚ₂ = float(first(l.cₚ)), float(last(l.cₚ))
            cₛ₁, cₛ₂ = float(first(l.cₛ)), float(last(l.cₛ))
            λₗ = max(cₚ₁, cₛ₁, cₚ₂, cₛ₂) / f
            # Kraken manual recommends double mesh density for elastic media
            ngl = round(Int, _value(2 * l.h / λₗ) * pm.mesh_density)
            push!(media, Medium{T}(depth, depth + l.h,
                [depth, depth + l.h],
                [crci(cₚ₁, in_dBperλ(l.δₚ), f; volume=fg),
                 crci(cₚ₂, in_dBperλ(l.δₚ), f; volume=fg)],
                [crci(cₛ₁, in_dBperλ(l.δₛ), f; volume=fg),
                 crci(cₛ₂, in_dBperλ(l.δₛ), f; volume=fg)],
                [ρ₁ / 1000, ρ₂ / 1000], false, ngl))
            push!(sigma, float(l.σ))
            depth += l.h
        end
        l = seabed.layers[end]
        seabed = ElasticBoundary(l.ρ, l.cₚ, l.cₛ, l.δₚ, l.δₛ)
    end

    hstop = _kj_halfspace(T, env.surface, f, fg)
    hsbot = _kj_halfspace(T, seabed, f, fg)
    push!(sigma, float(seabed === RigidBoundary || seabed === PressureReleaseBoundary ?
                       0.0 : seabed.σ))

    rmax = isinf(pm.rmax) ? T(1.01) * T(maxr) : T(pm.rmax)
    prob = Problem{T}(T(f), media, sigma, hstop, hsbot, T(pm.clow), T(pm.chigh),
                      rmax / 1000)      # rmax in km, as read from the env file

    # base mesh counts (NG); 0 → auto rule of ReadEnvironmentMod.f90
    ngs = [m.ng > 0 ? m.ng : auto_ng(m, f) for m in prob.media]
    prob, ngs
end

_kj_solve(pm::KrakenJL, tx, maxr) = begin
    prob, ngs = _make_problem(pm, tx, maxr)
    solve_modes(prob, ngs; complex_solver=pm.complex_solver, robust=pm.robust,
                threads=pm.threads)
end

### interface functions

function UnderwaterAcoustics.arrivals(pm::KrakenJL, tx1::AbstractAcousticSource,
                                      rx1::AbstractAcousticReceiver)
    res = _kj_solve(pm, tx1, location(rx1).x)
    ω = 2π * float(frequency(tx1))

    # group-speed validity check mirrors src/kraken.jl arrivals: KRAKEN (real
    # solver) leaves VG = 0 (reported as 0); invalid values → missing
    max_c = maximum(pm.env.soundspeed)
    if pm.env.seabed isa MultilayerElasticBoundary
        max_c = 2 * max(max_c, maximum(l -> max(maximum(l.cₚ), maximum(l.cₛ)),
                                       pm.env.seabed.layers[1:end-1]))
    elseif pm.env.seabed isa ElasticBoundary
        max_c = 2 * max(max_c, maximum(pm.env.seabed.cₚ), maximum(pm.env.seabed.cₛ))
    else
        max_c = 2 * max(max_c, maximum(pm.env.seabed.c))
    end
    T = real(eltype(res.k))
    v = all(vi -> 0 ≤ vi ≤ max_c, res.vg) ? Vector{Union{Missing,T}}(res.vg) :
        Vector{Union{Missing,T}}(fill(missing, length(res.k)))

    map(1:min(length(res.k), pm.nmodes)) do i
        # mode shapes on the solver's own mesh (finer than the wrapper's λ/10
        # resampling — see PORTING_NOTES §3)
        ψ = SampledField(res.phi[:, i]; z=-res.z)
        ModeArrival{Complex{T},typeof(ψ),Union{Missing,T},T}(
            i, res.k[i], ψ, v[i], ω / real(res.k[i]))
    end
end

function UnderwaterAcoustics.acoustic_field(pm::KrakenJL,
        tx1::AbstractAcousticSource, rx::AcousticReceiverGrid2D; mode=:coherent)
    mode ∈ (:coherent, :incoherent) || error("Unknown mode :" * string(mode))
    maxr = maximum(rx.xrange)
    res = _kj_solve(pm, tx1, maxr)
    M = min(length(res.k), pm.nmodes)

    zrev = first(rx.zrange) < last(rx.zrange)   # ascending depths for eval
    rxz = sort([-z for z in rx.zrange])
    xrev = first(rx.xrange) > last(rx.xrange)
    rxr = xrev ? reverse(collect(float.(rx.xrange))) : collect(float.(rx.xrange))

    phiS = tabulate_modes(res, [-location(tx1).z])[1, 1:M]
    phiR = tabulate_modes(res, rxz)[:, 1:M]
    P = evaluate_field(res.k[1:M], phiS, phiR, rxr;
                       incoherent=(mode === :incoherent), threads=pm.threads)

    # match the wrapper output: SHD files negate the pressure (_read_shd) and
    # Kraken uses the cis(-kᵣR) convention, hence the conjugation
    fld = conj.(-permutedims(P) .* db2amp(spl(tx1)))    # (nr, nz)
    xrev && (fld = reverse(fld; dims=1))
    zrev && (fld = reverse(fld; dims=2))
    fld
end

function UnderwaterAcoustics.acoustic_field(pm::KrakenJL,
        tx1::AbstractAcousticSource, rx1::AbstractAcousticReceiver; mode=:coherent)
    mode ∈ (:coherent, :incoherent) || error("Unknown mode :" * string(mode))
    p = location(rx1)
    res = _kj_solve(pm, tx1, p.x)
    M = min(length(res.k), pm.nmodes)
    phiS = tabulate_modes(res, [-location(tx1).z])[1, 1:M]
    phiR = tabulate_modes(res, [-p.z])[:, 1:M]
    P = evaluate_field(res.k[1:M], phiS, phiR, [float(p.x)];
                       incoherent=(mode === :incoherent), threads=pm.threads)
    conj(-P[1, 1] * db2amp(spl(tx1)))
end
