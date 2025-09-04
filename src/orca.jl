export Orca

"""
A propagation model based on the ORCA model.
"""
Base.@kwdef struct Orca{T} <: AbstractRayPropagationModel
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

function Orca(env; dz=0.1, kwargs...)
  _check_env(Orca, env)
  Orca{typeof(env)}(; env, kwargs...)
end

Base.show(io::IO, pm::Orca) = print(io, "Orca(⋯)")

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
      max(env.surface.c, 1e-6), max(env.surface.ρ / env.density, 1e-6), -in_dBperλ(env.surface.δ))
    ssp = env.soundspeed
    if is_constant(ssp)
      println(io, "*(3)\n2 0\n*(4)")
      @printf(io, "0 %0.6f 1 0\n", value(ssp))
      @printf(io, "%0.6f %0.6f\n", waterdepth, value(ssp))
    elseif ssp isa SampledFieldZ
      @printf(io, "*(3)\n%d 0\n*(4)\n", length(ssp.zrange))
      @printf(io, "%0.6f %0.6f 1 0\n", -z[1], ssp(z[1]))
      for z ∈ ssp.zrange[2:end]
        @printf(io, "%0.6f %0.6f\n", -z, ssp(z))
      end
    else
      n = _recommend_len(waterdepth, f)
      @printf(io, "*(3)\n%d 0\n*(4)\n", n)
      @printf(io, "%0.6f %0.6f 1 0\n", 0.0, ssp(0.0))
      for d ∈ range(0.0, waterdepth; length=n)[2:end]
        @printf(io, "%0.6f %0.6f\n", -z, ssp(z))
      end
    end
    if env.seabed isa MultilayerElasticBoundary
      @printf(io, "*(5)\n%d\n*(6)\n", length(env.seabed.layers)-1)
      for l ∈ env.seabed.layers[1:end-1]
        ρ₁, ρ₂ = first(l.ρ), last(l.ρ)
        cₚ₁, cₚ₂ = first(l.cₚ), last(l.cₚ)
        cₛ₁, cₛ₂ = first(l.cₛ), last(l.cₛ)
        aₚ = in_dBperλ(env.seabed.δₚ)
        aₛ = in_dBperλ(env.seabed.δₛ)
        @printf(io, "1 %0.6f %0.6f %0.6f %0.6f %0.6f %0.6f %0.6f %0.6f %0.6f %0.6f %0.6f 0 0 1 1\n",
          l.h, cₚ₁, cₚ₂, cₛ₁, cₛ₂, ρ₁, ρ₂, -aₚ, -aₚ, -aₛ, -aₛ)
        maxdepth += l.h
      end
      b = env.seabed.layers[end]
      @printf(io, "*(7)\n%0.6f %0.6f %0.6f %0.6f %0.6f 1 1\n",
        b.cₚ, b.cₛ, b.ρ / env.density, -in_dBperλ(b.δₚ), -in_dBperλ(b.δₛ))
    elseif env.seabed isa ElasticBoundary
      println(io, "*(5)\n0\n*(6)")
      @printf(io, "*(7)\n%0.6f %0.6f %0.6f %0.6f %0.6f 1 1\n",
        env.seabed.cₚ, env.seabed.cₛ, env.seabed.ρ / env.density,
        -in_dBperλ(env.seabed.δₚ), -in_dBperλ(env.seabed.δₛ))
    else
      println(io, "*(5)\n0\n*(6)")
      @printf(io, "*(7)\n%0.6f 0.0 %0.6f %0.6f 0.0 1 1\n",
        env.seabed.c, env.seabed.ρ / env.density, -in_dBperλ(env.seabed.δ))
    end
    println(io, "*(8)\n0\n*(9)")
  end
  open(opt_filename, "w") do io
    println(io, "*(1)\n2.01 1 0 0 0 0 3")
    @printf(io, "*(2)\n%d %0.6f %0.6f %0.6f %0.6f %0.6f %0.6f 0 0 0 0 0\n",
      pm.complex_solver ? 0 : 1, pm.cphmin, pm.cphmax, pm.rmin, pm.rmax, pm.phfac, pm.db_cut)
    @printf(io, "*(3)\n1 %0.6f\n", frequency(tx1))
    println(io, "*(4)\n1 1 0 0 0 1 0 0 0 0 0")
    @printf(io, "*(5)\n1 %0.6f\n", -location(tx1).z)
    if length(rx) == 1
      @printf(io, "1 %0.6f\n", -location(rx[1]).z)
      @printf(io, "1 %0.6f\n", location(rx[1]).x / 1000)
    elseif rx isa AcousticReceiverGrid2D
      x = rx.xrange ./ 1000
      z = -rx.zrange
      @printf(io, "%d %0.6f %0.6f\n", -length(z), minimum(z), maximum(z))
      @printf(io, "%d %0.6f %0.6f\n", -length(x), minimum(x), maximum(x))
    else
      error("Receivers must be on a 2D grid")
    end
    n = ceil(Int, maxdepth / pm.dz) + 1
    @printf(io, "*(6)\n0 1 %d 0 %0.6f\n", -n, maxdepth)
    println(io, "*(7)\n*(8)\n*(9)\n*(10)\n*(11)\n*(12)\n*(13)\n*(14)")
    range(0, maxdepth; length=n)
  end
end

### interface functions

function UnderwaterAcoustics.arrivals(pm::Orca, tx1::AbstractAcousticSource, rx1::AbstractAcousticReceiver)
  mktempdir(prefix="orca_") do dirname
    zs = _create_orca(pm, tx1, [rx1], dirname)
    _orca(dirname, pm.debug)
    ϕ = _read_modes_tlc(dirname).ϕ
    _read_orca_modes(dirname, zs, ϕ)
  end
end

function UnderwaterAcoustics.acoustic_field(pm::Orca, tx1::AbstractAcousticSource, rx::AcousticReceiverGrid2D)
  fld = mktempdir(prefix="orca_") do dirname
    _create_orca(pm, tx1, rx, dirname)
    _orca(dirname, pm.debug)
    fld = _read_modes_tlc(dirname).fld
  end
  fld .* db2amp(spl(tx1))
end

function UnderwaterAcoustics.acoustic_field(pm::Orca, tx1::AbstractAcousticSource, rx1::AbstractAcousticReceiver)
  fld = mktempdir(prefix="orca_") do dirname
    _create_orca(pm, tx1, [rx1], dirname)
    _orca(dirname, pm.debug)
    fld = _read_modes_tlc(dirname).fld
  end
  only(fld) .* db2amp(spl(tx1))
end

### helper functions

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
  ducts = Int[]
  modes = map(s[3:end]) do s1
    flds = split(strip(s1), r" +")
    i = parse(Int, flds[1])
    kre = parse(Float64, flds[2])
    att = parse(Float64, flds[3])
    kᵣ = ComplexF64(kre * Kw, log(10^(-att/1000/20)))
    vₚ = parse(Float64, flds[4])
    v = parse(Float64, flds[5])
    push!(ducts, parse(Int, flds[6]))
    ψ = SampledField(ϕ[:,i]; z=-zs)
    ModeArrival{ComplexF64,typeof(ψ),Union{Missing,Float64},Float64}(i, kᵣ, ψ, v, vₚ)
  end
  modes[ducts .== 1]
end

function _check_env(::Type{Orca}, env)
  env.seabed isa FluidBoundary || env.seabed isa ElasticBoundary || env.seabed isa MultilayerElasticBoundary || error("seabed must be a FluidBoundary, ElasticBoundary or MultilayerElasticBoundary")
  env.surface isa FluidBoundary || error("surface must be a FluidBoundary")
  is_range_dependent(env.soundspeed) && error("range-dependent soundspeed not supported")
  is_range_dependent(env.altimetry) && error("range-dependent altimetry not supported")
  is_range_dependent(env.bathymetry) && error("range-dependent bathymetry not supported")
  mktempdir(prefix="orca_") do dirname
    try
      _orca(dirname, false)
    catch e
      e isa ExecError && e.details == ["Unable to execute Orca"] && throw(e)
    end
  end
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
