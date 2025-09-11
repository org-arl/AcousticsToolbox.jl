export Bellhop

"""
A propagation model based on the FORTRAN OALIB Bellhop model.
"""
struct Bellhop{T} <: AbstractRayPropagationModel
  env::T
  nbeams::Int
  min_angle::Float32
  max_angle::Float32
  beam_type::Symbol
  beam_shift::Bool
  debug::Bool
  function Bellhop(env, nbeams, min_angle, max_angle, beam_type, beam_shift, debug)
    _check_env(Bellhop, env)
    -π/2 ≤ min_angle ≤ π/2 || error("min_angle should be between -π/2 and π/2")
    -π/2 ≤ max_angle ≤ π/2 || error("max_angle should be between -π/2 and π/2")
    min_angle < max_angle || error("max_angle should be more than min_angle")
    beam_type ∈ (:geometric, :gaussian) || error("Unknown beam_type type")
    new{typeof(env)}(env, max(nbeams, 0), Float32(in_units(u"rad", min_angle)), Float32(in_units(u"rad", max_angle)), beam_type, beam_shift, debug)
  end
end

"""
    Bellhop(env; kwargs...)

Create a Bellhop propagation model.

Supported keyword arguments:
- `nbeams`: number of beams to use (default: 0, auto)
- `min_angle`: minimum beam angle (default: -80°)
- `max_angle`: maximum beam angle (default: 80°)
- `beam_type`: `geometric` (default) or `gaussian`
- `beam_shift`: use beam shift (default: false)
- `debug`: debug mode (default: false)

Enabling debug mode will create a temporary directory with the Bellhop input and output files.
This allows manual inspection of the files.
"""
Bellhop(env; nbeams=0, min_angle=-80°, max_angle=80°, beam_type=:geometric, beam_shift=false, debug=false) = Bellhop(env, nbeams, min_angle, max_angle, beam_type, beam_shift, debug)

Base.show(io::IO, pm::Bellhop) = print(io, "Bellhop(⋯)")

### interface functions

function UnderwaterAcoustics.arrivals(pm::Bellhop, tx1::AbstractAcousticSource, rx1::AbstractAcousticReceiver; paths=true)
  mktempdir(prefix="bellhop_") do dirname
    nbeams = pm.nbeams
    nbeams == 0 && (nbeams = round(Int, (pm.max_angle - pm.min_angle) / deg2rad(0.05)) + 1)
    _write_env(pm, tx1, [rx1], dirname; nbeams, taskcode='A')
    _bellhop(dirname, pm.debug)
    arr = _read_arr(joinpath(dirname, "model.arr"))
    if paths
      _write_env(pm, tx1, [rx1], dirname; nbeams, taskcode='E')
      _bellhop(dirname, pm.debug)
      arr2 = _read_rays(joinpath(dirname, "model.ray"))
      for i ∈ eachindex(arr)
        ndx = findall(a -> a.ns == arr[i].ns && a.nb == arr[i].nb, arr2)
        if !isempty(ndx)
          ndx = argmin(j -> abs(arr2[j].θₛ - arr[i].θₛ), ndx)
          arr[i] = eltype(arr)(arr[i].t, arr[i].ϕ, arr[i].ns, arr[i].nb, arr[i].θₛ, arr[i].θᵣ, arr2[ndx].path)
        end
      end
    end
    arr
  end
end

function UnderwaterAcoustics.acoustic_field(pm::Bellhop, tx1::AbstractAcousticSource, rx::AcousticReceiverGrid2D; mode=:coherent)
  if mode === :coherent
    taskcode = 'C'
  elseif mode === :incoherent
    taskcode = 'I'
  elseif mode === :semicoherent
    taskcode = 'S'
  else
    error("Unknown mode :" * string(mode))
  end
  fld = mktempdir(prefix="bellhop_") do dirname
    xrev, zrev = _write_env(pm, tx1, rx, dirname; taskcode)
    _bellhop(dirname, pm.debug)
    _read_shd(joinpath(dirname, "model.shd"); xrev, zrev)
  end
  fld .* db2amp(spl(tx1))
end

function UnderwaterAcoustics.acoustic_field(pm::Bellhop, tx1::AbstractAcousticSource, rx1::AbstractAcousticReceiver; mode=:coherent)
  if mode === :coherent
    taskcode = 'C'
  elseif mode === :incoherent
    taskcode = 'I'
  elseif mode === :semicoherent
    taskcode = 'S'
  else
    error("Unknown mode :" * string(mode))
  end
  fld = mktempdir(prefix="bellhop_") do dirname
    _write_env(pm, tx1, [rx1], dirname; taskcode)
    _bellhop(dirname, pm.debug)
    _read_shd(joinpath(dirname, "model.shd"))[1]
  end
  fld .* db2amp(spl(tx1))
end

"""
    rays(pm::Bellhop, tx, max_range)
    rays(pm::Bellhop, tx, max_range, nbeams)

Compute ray trace for a single transmitter.
"""
function rays(pm::Bellhop, tx1::AbstractAcousticSource, max_range::Real, nbeams::Int=pm.nbeams)
  mktempdir(prefix="bellhop_") do dirname
    _write_env(pm, tx1, [AcousticReceiver(max_range, 0)], dirname; nbeams, taskcode='R')
    _bellhop(dirname, pm.debug)
    _read_rays(joinpath(dirname, "model.ray"))
  end
end

### helpers

function _bellhop(dirname, debug)
  infilebase = joinpath(dirname, "model")
  outfilename = joinpath(dirname, "output.txt")
  try
    run(pipeline(ignorestatus(`$(AcousticsToolbox_jll.bellhop()) $infilebase`); stdout=outfilename, stderr=outfilename))
    if debug
      @info "Bellhop run completed in $dirname, press ENTER to delete intermediate files..."
      readline()
    end
  catch
    throw(ExecError("Bellhop", ["Unable to execute Bellhop"]))
  end
  err = String[]
  _check_err!(err, outfilename)
  _check_err!(err, joinpath(dirname, "model.prt"))
  if length(err) > 0
    throw(ExecError("Bellhop", err))
  end
end

function _check_env(::Type{Bellhop}, env)
  env.seabed isa FluidBoundary || env.seabed isa ElasticBoundary || error("seabed must be a FluidBoundary or ElasticBoundary")
  env.surface isa FluidBoundary || error("surface must be a FluidBoundary")
  is_range_dependent(env.soundspeed) && error("Range-dependent soundspeed not supported")
  nothing
end
