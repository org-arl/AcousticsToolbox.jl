export Orca

"""
A propagation model based on the ORCA model.
"""
Base.@kwdef struct Orca{T} <: AbstractModePropagationModel
  env::T
  dz::Float64 = 0.1
  complex_solver::Bool = true
  cphmin::Float64 = 0.0
  cphmax::Float64 = 0.0
  rmin::Float64 = 0.0
  rmax::Float64 = 0.0
  phfac::Float64 = 0.0
  db_cut::Float64 = 0.0
  debug::Bool = false
end

"""
    Orca(env; kwargs...)

Create a Orca propagation model.

Supported keyword arguments:
- `dz`: vertical step size for sampling modes in m (default: 0.1)
- `complex_solver`: selects between real/complex solver (default: `true`)
- `cphmin`: minimum phase speed in m/s (default: 0, 0=auto)
- `cphmax`: maximum phase speed in m/s (default: 0, 0=auto)
- `rmin`: minimum range of interest in m (default: 0, 0=auto)
- `rmax`: maximum range of interest in m (default: 0, 0=auto)
- `phfac`: phase step parameter (default: 0, 0=auto)
- `db_cut`: dB cutoff for weak modes (default: 0, 0=auto)
- `debug`: debug mode (default: false)

Enabling debug mode will create a temporary directory with the Orca input and output files.
This allows manual inspection of the files.
"""
function Orca(env; kwargs...)
  _check_env(Orca, env)
  Orca{typeof(env)}(; env, kwargs...)
end

Base.show(io::IO, pm::Orca) = print(io, "Orca(⋯)")

### interface functions

function UnderwaterAcoustics.arrivals(pm::Orca, tx1::AbstractAcousticSource, rx1::AbstractAcousticReceiver)
  mktempdir(prefix="orca_") do dirname
    D, n = _create_orca(pm, tx1, [rx1], dirname)::Tuple{Float64,Int}
    zs = range(0, D; length=n)
    _orca(dirname, pm.debug)
    ϕ = _read_modes_tlc(dirname).ϕ
    _read_orca_modes(dirname, zs, ϕ)
  end
end

function UnderwaterAcoustics.acoustic_field(pm::Orca, tx1::AbstractAcousticSource, rx::AcousticReceiverGrid2D)
  fld = mktempdir(prefix="orca_") do dirname
    _create_orca(pm, tx1, rx, dirname)
    _orca(dirname, pm.debug)
    x = _read_modes_tlc(dirname).fld
    # HACK: extra run because the first TL in the array from ORCA is incorrect
    if size(rx, 1) > 1
      rx1 = AcousticReceiverGrid2D(rx.xrange[1], rx.zrange)
      _create_orca(pm, tx1, rx1, dirname)
      _orca(dirname, pm.debug)
      y = _read_modes_tlc(dirname).fld
      x[1,:] .= y[2,:]
    else
      x = x[2:end,:]
    end
    x
  end
  fld .* db2amp(spl(tx1))
end

function UnderwaterAcoustics.acoustic_field(pm::Orca, tx1::AbstractAcousticSource, rx1::AbstractAcousticReceiver)
  fld = mktempdir(prefix="orca_") do dirname
    _create_orca(pm, tx1, [rx1], dirname)
    _orca(dirname, pm.debug)
    _read_modes_tlc(dirname).fld
  end
  fld[2] .* db2amp(spl(tx1))
end

### helper functions

function _create_orca(pm, tx1, rx, dirname)
  location(tx1).x == 0 || error("Transmitter must be at (0, 0, z)")
  location(tx1).y == 0 || error("2D model requires transmitter in the x-z plane")
  all(location(rx1).x >= 0 for rx1 ∈ rx) || error("Receivers must be in the +x halfspace")
  all(location(rx1).y == 0 for rx1 ∈ rx) || error("2D model requires receivers in the x-z plane")
  name = basename(dirname)
  svp_filename = joinpath(dirname, "orca.svp")
  opt_filename = joinpath(dirname, "orca.opt")
  orca_filename = joinpath(dirname, "orca_in")
  open(orca_filename, "w") do io
    println(io, basename(svp_filename))
    println(io, basename(opt_filename))
    println(io, "orca.out")
  end
  env = pm.env
  waterdepth = maximum(env.bathymetry)
  maxdepth = waterdepth
  open(svp_filename, "w") do io
    println(io, "*(1)\n3.0 '$name'")
    @printf(io, "*(2)\n%0.6f 0.0 %0.6f %0.6f 0.0 1.0 1.0\n",
      _c(env.surface.c), _ρ(env.surface.ρ / env.density), -in_dBperλ(env.surface.δ))
    ssp = env.soundspeed
    f = frequency(tx1)
    att = -20*log10(absorption(f, 1, pm.env.salinity, pm.env.temperature, waterdepth/2, pm.env.pH)) / f * 1000
    if is_constant(ssp)
      println(io, "*(3)\n2 0\n*(4)")
      @printf(io, "0 %0.6f 1 %0.6f\n", value(ssp), att)
      @printf(io, "%0.6f %0.6f\n", waterdepth, value(ssp))
    elseif ssp isa SampledFieldZ
      @printf(io, "*(3)\n%d 0\n*(4)\n", length(ssp.zrange))
      @printf(io, "%0.6f %0.6f 1 %0.6f\n", -ssp.zrange[1], ssp(ssp.zrange[1]), att)
      for z ∈ ssp.zrange[2:end]
        @printf(io, "%0.6f %0.6f\n", -z, ssp(z))
      end
    else
      n = _recommend_len(waterdepth, f)
      @printf(io, "*(3)\n%d 0\n*(4)\n", n)
      @printf(io, "%0.6f %0.6f 1 %0.6f\n", 0.0, ssp(0.0), att)
      for d ∈ range(0.0, waterdepth; length=n)[2:end]
        @printf(io, "%0.6f %0.6f\n", -d, ssp(d))
      end
    end
    if env.seabed isa MultilayerElasticBoundary
      @printf(io, "*(5)\n%d\n*(6)\n", length(env.seabed.layers)-1)
      for l ∈ env.seabed.layers[1:end-1]
        ρ₁, ρ₂ = first(l.ρ), last(l.ρ)
        cₚ₁, cₚ₂ = first(l.cₚ), last(l.cₚ)
        cₛ₁, cₛ₂ = first(l.cₛ), last(l.cₛ)
        aₚ = in_dBperλ(l.δₚ)
        aₛ = in_dBperλ(l.δₛ)
        @printf(io, "1 %0.6f %0.6f %0.6f %0.6f %0.6f %0.6f %0.6f %0.6f %0.6f %0.6f %0.6f 0 0 1 1\n",
          l.h, _c(cₚ₁), _c(cₚ₂), cₛ₁, cₛ₂, _ρ(ρ₁ / env.density), _ρ(ρ₂ / env.density), -aₚ, -aₚ, -aₛ, -aₛ)
        maxdepth += l.h
      end
      b = env.seabed.layers[end]
      @printf(io, "*(7)\n%0.6f %0.6f %0.6f %0.6f %0.6f 1 1\n",
        _c(b.cₚ), b.cₛ, _ρ(b.ρ / env.density), -in_dBperλ(b.δₚ), -in_dBperλ(b.δₛ))
    elseif env.seabed isa ElasticBoundary
      println(io, "*(5)\n0\n*(6)")
      @printf(io, "*(7)\n%0.6f %0.6f %0.6f %0.6f %0.6f 1 1\n",
        _c(env.seabed.cₚ), env.seabed.cₛ, _ρ(env.seabed.ρ / env.density),
        -in_dBperλ(env.seabed.δₚ), -in_dBperλ(env.seabed.δₛ))
    else
      println(io, "*(5)\n0\n*(6)")
      @printf(io, "*(7)\n%0.6f 0.0 %0.6f %0.6f 0.0 1 1\n",
        _c(env.seabed.c), _ρ(env.seabed.ρ / env.density), -in_dBperλ(env.seabed.δ))
    end
    println(io, "*(8)\n0\n*(9)")
  end
  open(opt_filename, "w") do io
    println(io, "*(1)\n2.01 1 0 0 0 0 3")
    @printf(io, "*(2)\n%d %0.6f %0.6f %0.6f %0.6f %0.6f %0.6f 0 0 0 0 0\n",
      pm.complex_solver ? 0 : 1, pm.cphmin, pm.cphmax, pm.rmin / 1000, pm.rmax / 1000, pm.phfac, pm.db_cut)
    @printf(io, "*(3)\n1 %0.6f\n", frequency(tx1))
    println(io, "*(4)\n1 1 0 0 0 1 0 0 0 0 0")
    @printf(io, "*(5)\n1 %0.6f\n", -location(tx1).z)
    if length(rx) == 1
      @printf(io, "1 %0.6f\n", -location(rx[1]).z)
      # HACK: extra entry at start because the first TL from ORCA is incorrect
      x = location(rx[1]).x / 1000
      @printf(io, "2 %0.6f %0.6f\n", x, x)
    elseif rx isa AcousticReceiverGrid2D
      x = rx.xrange ./ 1000
      z = -rx.zrange
      @printf(io, "%d %0.6f %0.6f\n", -length(z), minimum(z), maximum(z))
      if size(rx, 1) == 1
        # HACK: extra entry at start because the first TL from ORCA is incorrect
        @printf(io, "2 %0.6f %0.6f\n", x[1,1], x[1,1])
      else
        @printf(io, "%d %0.6f %0.6f\n", -length(x), minimum(x), maximum(x))
      end
    else
      error("Receivers must be on a 2D grid")
    end
    n = min(ceil(Int, maxdepth / pm.dz) + 1, 10000)::Int64
    @printf(io, "*(6)\n0 1 %d 0 %0.6f\n", -n, maxdepth)
    println(io, "*(7)\n*(8)\n*(9)\n*(10)\n*(11)\n*(12)\n*(13)\n*(14)")
    (Float64(maxdepth), n)
  end
end

_c(c) = clamp(c, 1e-6, 1e8)
_ρ(ρ) = clamp(ρ, 1e-6, 1e8)

function _orca(dirname, debug)
  outfilename = joinpath(dirname, "output.txt")
  try
    run(pipeline(Cmd(AcousticsToolbox_jll.orca(); ignorestatus=true, dir=dirname); stdout=outfilename, stderr=outfilename))
    if debug
      @info "Orca run completed in $dirname, press ENTER to delete intermediate files..."
      readline()
    end
  catch
    throw(ExecError("Orca", ["Unable to execute Orca"]))
  end
end

function _read_orca_modes(dirname, zs, ϕ)
  filename = joinpath(dirname, "orca.out_modes")
  s = readlines(filename)
  Kw = parse(Float64, match(r"Kw = ([\d\.E\+\-]+)", s[1])[1])
  map(3:length(s)) do i
    flds = split(strip(s[i]), r" +")
    m = parse(Int, flds[1])
    kre = parse(Float64, flds[2])
    att = parse(Float64, flds[3])
    kᵣ = ComplexF64(kre * Kw, log(10^(-att/1000/20)))
    vₚ = parse(Float64, flds[4])
    v = parse(Float64, flds[5])
    ψ = SampledField(ϕ[:,m]; z=-zs)
    ModeArrival{ComplexF64,typeof(ψ),Float64,Float64}(m, kᵣ, ψ, v, vₚ)
  end
end

function _check_env(::Type{Orca}, env)
  env.seabed isa FluidBoundary || env.seabed isa ElasticBoundary || env.seabed isa MultilayerElasticBoundary || error("seabed must be a FluidBoundary, ElasticBoundary or MultilayerElasticBoundary")
  env.surface isa FluidBoundary || error("surface must be a FluidBoundary")
  is_range_dependent(env.soundspeed) && error("Range-dependent soundspeed not supported")
  is_range_dependent(env.altimetry) && error("Range-dependent altimetry not supported")
  is_range_dependent(env.bathymetry) && error("Range-dependent bathymetry not supported")
  nothing
end

function _read_modes_tlc(dirname)
  filename = joinpath(dirname, "modes_tlc.bin")
  open(filename, "r") do io
    tag = read(io, UInt32)
    nzsr = read(io, UInt32)
    nmode = read(io, UInt32)
    nrec = read(io, UInt32)
    nsrc = read(io, UInt32)
    read(io, UInt32) == tag || error("Bad modes output")
    tag = read(io, UInt32)
    phi = Array{ComplexF32}(undef, nzsr, nmode)
    for jm in 1:nmode
      for j in 1:nzsr
        phi[j, jm] = read(io, ComplexF32)
      end
    end
    read(io, UInt32) == tag || error("Bad modes output")
    tag = read(io, UInt32)
    tlc = Array{ComplexF32}(undef, nrec, nsrc)
    for jsrc in 1:nsrc
      for jrec in 1:nrec
        tlc[jrec, jsrc] = read(io, ComplexF32)
      end
    end
    read(io, UInt32) == tag || error("Bad modes output")
    (ϕ=phi, fld=reverse(transpose(tlc); dims=2))
  end
end
