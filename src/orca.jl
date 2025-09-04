using SimpleDiffEq

export Orca

"""
A propagation model based on the ORCA model.
"""
Base.@kwdef struct Orca{T} <: AbstractRayPropagationModel
  env::T
  complex_solver::Bool = true
  clow::Float32 = 0.0
  chigh::Float32 = 0.0
  rmin::Float32 = 0.0
  rmax::Float32 = 0.0
  phfac::Float32 = 1.0
  cutoff::Float32 = 0.0
  debug::Bool = false
end

Orca(env; kwargs...) = Orca{typeof(env)}(; env, kwargs...)

Base.show(io::IO, pm::Orca) = print(io, "Orca(⋯)")

function _create_orca(pm, tx, rx, dirname)
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
  open(svp_filename, "w") do io
    println(io, "*(1)\n2.0 '$name'")
    @printf(io, "*(2)\n%0.6f 0.0 %0.6f %0.6f 0.0\n", max(env.surface.c, 1e-6), max(env.surface.ρ / env.density, 1e-6), -in_dBperλ(env.surface.δ))
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
        @printf(io, "1 %0.6f %0.6f %0.6f %0.6f %0.6f %0.6f %0.6f %0.6f %0.6f %0.6f\n",
          cₚ₁, cₚ₂, cₛ₁, cₛ₂, ρ₁, ρ₂, -env.seabed.δₚ, -env.seabed.δₚ, -env.seabed.δₛ, -env.seabed.δₛ)
      end
      b = env.seabed.layers[end]
      @printf(io, "*(7)\n%0.6f %0.6f %0.6f %0.6f %0.6f\n", b.cₚ, b.cₛ, b.ρ / env.density, -b.δₚ, -b.δₛ)
    elseif env.seabed isa ElasticBoundary
      println(io, "*(5)\n0\n*(6)")
      @printf(io, "*(7)\n%0.6f %0.6f %0.6f %0.6f %0.6f\n", env.seabed.cₚ, env.seabed.cₛ, env.seabed.ρ / env.density, -env.seabed.δₚ, -env.seabed.δₛ)
    else
      println(io, "*(5)\n0\n*(6)")
      @printf(io, "*(7)\n%0.6f 0.0 %0.6f %0.6f 0.0\n", env.seabed.c, env.seabed.ρ / env.density, -env.seabed.δ)  # TODO check att units
    end
    println(io, "*(8)\n0\n*(9)")
  end
  open(opt_filename, "w") do io
    println(io, "*(1)\n1.6 1 0 0 0 0 3")
    @printf(io, "*(2)\n%d %0.6f %0.6f %0.6f %0.6f %0.6f %0.6f 0 0 0\n",
      pm.complex_solver ? 0 : 1, pm.clow, pm.chigh, pm.rmin, pm.rmax, pm.phfac, pm.cutoff)
    flist = [frequency(tx1) for tx1 ∈ tx]
    f = sum(flist) / length(flist)
    maximum(abs, flist .- f) / f > 0.2 && @warn("Source frequency varies by more than 20% from nominal frequency")
    @printf(io, "*(3)\n1 %0.6f\n", f)
    println(io, "*(4)\n1 0 0 0 3 1 0 0 0 0 0")
    @printf(io, "*(5)\n%d", length(tx))
    for tx1 ∈ tx
      @printf(io, " %0.6f", -location(tx1).z)
    end
    if length(rx) == 1
      @printf(io, "\n1 %0.6f\n", -location(rx[1]).z)
      @printf(io, "1 %0.6f\n", location(rx[1]).x / 1000)
    elseif rx isa AcousticReceiverGrid2D
      x = rx.xrange ./ 1000
      z = -rx.zrange
      @printf(io, "\n%d %0.6f %0.6f\n", -length(z), minimum(z), maximum(z))
      @printf(io, "%d %0.6f %0.6f\n", -length(x), minimum(x), maximum(x))
    else
      error("Receivers must be on a 2D grid")
    end
    println(io, "*(6)\n*(7)\n*(8)\n*(9)\n*(10)\n*(11)\n*(12)\n*(13)\n*(14)")
  end
end

### interface functions

function UnderwaterAcoustics.arrivals(pm::Orca, tx1::AbstractAcousticSource, rx1::AbstractAcousticReceiver)
  mktempdir(prefix="orca_") do dirname
    _create_orca(pm, [tx1], [rx1], dirname)
    _orca(dirname, pm.debug)
    _read_orca_modes(dirname, frequency(tx1), pm.env.soundspeed, pm.env.bathymetry, 0.1) #pm.dz) # FIXME: extend depth to sediment layers
  end
end

function UnderwaterAcoustics.acoustic_field(pm::Orca, tx1::AbstractAcousticSource, rx::AcousticReceiverGrid2D)
  fld = mktempdir(prefix="orca_") do dirname
    _create_orca(pm, [tx1], rx, dirname)
    _orca(dirname, pm.debug)
    # TODO
  end
  # fld .* db2amp(spl(tx1))
end

function UnderwaterAcoustics.acoustic_field(pm::Orca, tx1::AbstractAcousticSource, rx1::AbstractAcousticReceiver)
  fld = mktempdir(prefix="orca_") do dirname
    _create_orca(pm, [tx1], [rx1], dirname)
    _orca(dirname, pm.debug)
    # TODO
  end
  # fld .* db2amp(spl(tx1))
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

function _read_orca_modes(dirname, f, c, D, dz)
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
    ψ = _compute_mode_function(f, c, kᵣ, D, dz)
    ModeArrival{ComplexF64,typeof(ψ),Union{Missing,Float64},Float64}(i, kᵣ, ψ, v, vₚ)
  end
  modes[ducts .== 1]
end

# FIXME: assumes pressure release surface and no density variation
function _compute_mode_function(f, c, kᵣ, D, dz)
  ω = 2π * f
  function du!(du, u, p, t)
    du[1] = u[2]
    du[2] = -((ω / value(c, -t))^2 - kᵣ^2) * u[1]
  end
  prob = ODEProblem{true}(du!, ComplexF64[0, 1], (0, D))
  sol = solve(prob, SimpleTsit5(); dt=dz)
  zs = range(0, -D; step=-dz)
  ψs = map(z -> first(sol(-z)), zs)
  ψs ./= sqrt(sum(abs2, ψs) * dz)
  SampledField(ψs; z=zs)
end
