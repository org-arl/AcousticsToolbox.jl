export Kraken

const KRAKEN = Ref{Cmd}(AcousticsToolbox_jll.kraken())
const KRAKENC = Ref{Cmd}(AcousticsToolbox_jll.krakenc())
const FIELD = Ref{Cmd}(AcousticsToolbox_jll.field())

"""
A propagation model based on the FORTRAN OALIB Kraken model.
"""
struct Kraken{T} <: AbstractModePropagationModel
  env::T
  nmodes::Int
  nmesh::Int
  clow::Float32
  chigh::Float32
  leaky::Bool
  debug::Bool
  function Kraken(env, nmodes, nmesh, clow, chigh, leaky, debug)
    _check_env(Kraken, env)
    nmodes ≥ 1 || throw(ArgumentError("number of modes should be positive"))
    nmesh ≥ 0 || throw(ArgumentError("number of mesh points should be non-negative"))
    clow ≥ 0.0 || throw(ArgumentError("clow should be non-negative"))
    chigh > clow || throw(ArgumentError("chigh should be more than clow"))
    new{typeof(env)}(env, nmodes, nmesh, clow, chigh, leaky, debug)
  end
end

"""
    Kraken(env; kwargs...)

Create a Kraken propagation model.

Supported keyword arguments:
- `nmodes`: number of modes to use (default: 9999)
- `nmesh`: number of mesh points (default: 0, auto)
- `clow`: lower limit of phase speed (default: 0, auto)
- `chigh`: upper limit of phase speed (default: 1600.0)
- `leaky`: use KrakenC for leaky modes (default: false)
- `debug`: debug mode (default: false)

Enabling debug mode will create a temporary directory with the Kraken input and output files.
This allows manual inspection of the files.
"""
Kraken(env; nmodes=9999, nmesh=0, clow=0.0, chigh=1600.0, leaky=false, debug=false) = Kraken(env, nmodes, nmesh, clow, chigh, leaky, debug)

Base.show(io::IO, pm::Kraken) = print(io, "Kraken(⋯)")

### interface functions

# function UnderwaterAcoustics.arrivals(pm::Kraken, tx1::AbstractAcousticSource, rx1::AbstractAcousticReceiver)
#   mktempdir(prefix="kraken_") do dirname
#     # TODO
#   end
# end

function UnderwaterAcoustics.acoustic_field(pm::Kraken, tx1::AbstractAcousticSource, rx::AcousticReceiverGrid2D; mode=:coherent)
  if mode === :coherent
    taskcode = 'C'
  elseif mode === :incoherent
    taskcode = 'I'
  else
    throw(ArgumentError("Unknown mode :" * string(mode)))
  end
  mktempdir(prefix="kraken_") do dirname
    xrev, zrev = _write_env(pm, [tx1], rx, dirname)
    _kraken(dirname, pm.leaky, pm.debug)
    _write_flp(pm, [tx1], rx, dirname)
    _field(dirname, pm.debug)
    _read_shd(joinpath(dirname, "model.shd"); xrev, zrev)
  end
end

function UnderwaterAcoustics.acoustic_field(pm::Bellhop, tx1::AbstractAcousticSource, rx1::AbstractAcousticReceiver; mode=:coherent)
  if mode === :coherent
    taskcode = 'C'
  elseif mode === :incoherent
    taskcode = 'I'
  else
    throw(ArgumentError("Unknown mode :" * string(mode)))
  end
  mktempdir(prefix="bellhop_") do dirname
    _write_env(pm, [tx1], [rx1], dirname)
    _kraken(dirname, pm.leaky, pm.debug)
    _write_flp(pm, [tx1], [rx1], dirname)
    _field(dirname, pm.debug)
    _read_shd(joinpath(dirname, "model.shd"))[1]
  end
end

### helper functions

function _check_env(::Type{Kraken}, env)
  env.seabed isa FluidBoundary || throw(ArgumentError("seabed must be a FluidBoundary"))
  env.surface isa FluidBoundary || throw(ArgumentError("surface must be a FluidBoundary"))
  is_range_dependent(env.soundspeed) && throw(ArgumentError("range-dependent soundspeed not supported"))
  mktempdir(prefix="kraken_") do dirname
    try
      kraken(dirname, false)
      krakenc(dirname, false)
    catch e
      e isa ExecError && e.details == ["Unable to execute Kraken"] && throw(e)
    end
  end
  nothing
end

function _kraken(dirname, leaky, debug)
  infilebase = joinpath(dirname, "model")
  outfilename = joinpath(dirname, "output.txt")
  try
    run(pipeline(ignorestatus(`$(leaky ? KRAKENC[] : KRAKEN[]) $infilebase`); stdout=outfilename, stderr=outfilename))
    if debug
      @info "Kraken run completed in $dirname, press ENTER to delete intermediate files..."
      readline()
    end
  catch
    throw(ExecError("Kraken", ["Unable to execute Kraken"]))
  end
  err = String[]
  _check_err!(err, outfilename)
  _check_err!(err, joinpath(dirname, "model.prt"))
  if length(err) > 0
    throw(ExecError("Kraken", err))
  end
end

function _field(dirname, debug)
  infilebase = joinpath(dirname, "model")
  outfilename = joinpath(dirname, "output.txt")
  try
    run(pipeline(ignorestatus(Cmd(`$(FIELD[]) $infilebase`; dir=dirname)); stdout=outfilename, stderr=outfilename))
    if debug
      @info "Field run completed in $dirname, press ENTER to delete intermediate files..."
      readline()
    end
  catch
    throw(ExecError("Field", ["Unable to execute Field"]))
  end
  err = String[]
  _check_err!(err, outfilename)
  _check_err!(err, joinpath(dirname, "field.prt"))
  if length(err) > 0
    throw(ExecError("Field", err))
  end
end

function _write_flp(pm::Kraken, tx::Vector{<:AcousticSource}, rx::AbstractArray{<:AcousticReceiver}, dirname)
  filename = joinpath(dirname, "model.flp")
  name = split(basename(dirname), "_")[end]
  open(filename, "w") do io
    println(io, "'", name, "'")        # 1) title
    length(tx) > 1 && (src = "X")      # 2) line source
    length(tx) == 1 && (src = "R")     # 2) point source
    # 2)  Adiabatic mode thoery TODO: add coupled model
    # TODO: check why coherent option cause error
    print(io, "'", src, "A", "'\n")    # 2) Coherent / incoherent mode addition
    @printf(io, "%i\n", pm.nmodes)     # 3) Number of modes
    println(io, "1")                   # 4) Number of profile    TODO: hardcoded, check for general setting
    println(io, "0.0")                 # 4) Profile ranges (km)  TODO: hardcoded, check for general setting
    if length(rx) == 1
      _print_array(io, [location(rx[1])[1]./ 1000.0] )     # 6) receiver ranges(km)
      _print_array(io, [-location(tx1)[3] for tx1 ∈ tx])   #    source depth (m) TODO: check multiple sources how
      _print_array(io, [-location(rx1)[3] for rx1 ∈ rx] )  # 6) receiver depth(m)
      nr = 1
    elseif rx isa AcousticReceiverGrid2D
      d = reverse(-rx.zrange)
      if first(d) > last(d)
       d = reverse(d)
       zrev = true
      end
      r = rx.xrange ./ 1000.0
      if first(r) > last(r)
       r = reverse(r)
       xrev = true
      end
      _print_array(io, r)                                 #    receiver ranges (km)
      _print_array(io, [-location(tx1)[3] for tx1 ∈ tx])  #    source depth (m)
      _print_array(io, d)                                 # 6) receiver depth(m)
      nr = length(d)
    end
    _print_array(io, zeros(nr))                           # number of receiver range displacement (same as NRD)
  end
end
