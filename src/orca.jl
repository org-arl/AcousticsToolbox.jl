export Orca

const ORCA = Ref{Cmd}(`orca90.exe`) # (AcousticsToolbox_jll.orca90())

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
  false_bottom::Bool = false
  mode_depths::Int = 0
end

function _create_orca(pm, tx, rx, dirname)
  name = basename(dirname)
  svp_filename = joinpath(dirname, "orca.svp")
  opt_filename = joinpath(dirname, "orca.opt")
  out_filename = joinpath(dirname, "orca.out")
  orca_filename = joinpath(dirname, "orca_in")
  open(orca_filename, "w") do io
    println(io, svp_filename)
    println(io, opt_filename)
    println(io, out_filename)
  end
  env = pm.env
  waterdepth = maximum(env.bathymetry)
  open(svp_filename, "w") do io
    println(io, "*\n3.0 '$name'")
    @printf(io, "*\n%0.6f 0.0 %0.6f %0.6f 0.0 1.0 1.0\n", env.surface.c, env.surface.ρ / env.density, -in_dBperλ(env.surface.δ))
    ssp = env.soundspeed
    if is_constant(ssp)
      println(io, "*\n2 0\n*")
      @printf(io, "0 %0.6f 1 999\n", value(ssp))
      @printf(io, "%0.6f %0.6f\n", waterdepth, value(ssp))
    elseif ssp isa SampledFieldZ
      @printf(io, "*\n%d 0\n*\n", length(ssp.zrange))
      @printf(io, "%0.6f %0.6f 1 999\n", -z[1], ssp(z[1]))
      for z ∈ ssp.zrange[2:end]
        @printf(io, "%0.6f %0.6f\n", -z, ssp(z))
      end
    else
      n = _recommend_len(waterdepth, f)
      @printf(io, "*\n%d 0\n*\n", n)
      @printf(io, "%0.6f %0.6f 1 999\n", 0.0, ssp(0.0))
      for d ∈ range(0.0, waterdepth; length=n)[2:end]
        @printf(io, "%0.6f %0.6f\n", -z, ssp(z))
      end
    end
    println(io, "*\n0\n*")   # TODO: support bottom layers
    @printf(io, "*\n%0.6f 0.0 %0.6f %0.6f 0.0 1.0 1.0\n", env.seabed.c, env.seabed.ρ / env.density, -in_dBperλ(env.seabed.δ))   # TODO: support elastic seabed
    println(io, "*\n0\n*")
  end
  open(opt_filename, "w") do io
    println(io, "*\n1.6 1 0 0 0 0 0")   # TODO: change last param to choose specific output format
    @printf(io, "*\n%d %0.6f %0.6f %0.6f %0.6f %0.6f %0.6f %d 0 0\n",
      pm.leaky ? 0 : 1, pm.min_phase_speed, pm.max_phase_speed, pm.rmin,
      pm.rmax, pm.phase_step, pm.cutoff, pm.false_bottom ? 1 : 0)
    flist = [frequency(tx1) for tx1 ∈ tx]
    f = sum(flist) / length(flist)
    maximum(abs, flist .- f) / f > 0.2 && @warn("Source frequency varies by more than 20% from nominal frequency")
    @printf(io, "*\n1 %0.6f\n", f)
    @printf(io, "*\n1 %d 0 0 0 1 1 0 0 0 0\n", pm.mode_depths > 0 ? 1 : 0)
    @printf(io, "*\n%d", length(tx))
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
      @printf(io, "*\n3 0 %d 0 %0.6f\n", -pm.mode_depths, waterdepth)
    else
      println(io, "*")
    end
    println(io, "*\n*\n*\n*\n*\n*\n*\n*")
  end
end
