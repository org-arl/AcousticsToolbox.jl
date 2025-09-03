export Orca

"""
A propagation model based on the ORCA model.
"""
Base.@kwdef struct Orca{T} <: AbstractRayPropagationModel
  env::T
  leaky::Bool = true
  min_phase_speed::Float32 = 0.0
  max_phase_speed::Float32 = 0.0
  rmin::Float32 = 0.0
  rmax::Float32 = 0.0
  phase_step::Float32 = 1.0
  cutoff::Float32 = 0.0
  false_bottom::Bool = false          # TODO: remove
  mode_depths::Int = 0
  debug::Bool = false
end

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
    @printf(io, "*(2)\n%0.6f 0.0 %0.6f %0.6f 0.0\n", env.surface.c, env.surface.ρ / env.density, -in_dBperλ(env.surface.δ))
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
    println(io, "*(5)\n0\n*(6)")   # TODO: support bottom layers
    @printf(io, "*(7)\n%0.6f 0.0 %0.6f %0.6f 0.0\n", env.seabed.c, env.seabed.ρ / env.density, -in_dBperλ(env.seabed.δ))   # TODO: support elastic seabed
    println(io, "*(8)\n0\n*(9)")
  end
  open(opt_filename, "w") do io
    println(io, "*(1)\n2.01 1 0 0 0 0 3")
    @printf(io, "*(2)\n%d %0.6f %0.6f %0.6f %0.6f %0.6f %0.6f 0 0 0 0 0\n",
      pm.leaky ? 0 : 1, pm.min_phase_speed, pm.max_phase_speed, pm.rmin,
      pm.rmax, pm.phase_step, pm.cutoff)
    flist = [frequency(tx1) for tx1 ∈ tx]
    f = sum(flist) / length(flist)
    maximum(abs, flist .- f) / f > 0.2 && @warn("Source frequency varies by more than 20% from nominal frequency")
    @printf(io, "*(3)\n1 %0.6f\n", f)
    @printf(io, "*(4)\n1 %d 0 0 0 1 1 0 0 0 0\n", pm.mode_depths > 0 ? 1 : 0)
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
    if pm.mode_depths > 0
      @printf(io, "*(6)\n3 0 %d 0 %0.6f\n", -pm.mode_depths, waterdepth)
    else
      println(io, "*(6)")
    end
    println(io, "*(7)\n*(8)\n*(9)\n*(10)\n*(11)\n*(12)\n*(13)\n*(14)")
  end
end

### interface functions

function UnderwaterAcoustics.arrivals(pm::Orca, tx1::AbstractAcousticSource, rx1::AbstractAcousticReceiver)
  mktempdir(prefix="orca_") do dirname
    _create_orca(pm, [tx1], [rx1], dirname)
    _orca(dirname, pm.debug)
    # TODO
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
