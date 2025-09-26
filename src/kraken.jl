export Kraken

"""
A propagation model based on the FORTRAN OALIB Kraken model.
"""
struct Kraken{T} <: AbstractModePropagationModel
  env::T
  nmodes::Int
  mesh_density::Float32
  clow::Float32
  chigh::Float32
  rmax::Float32
  complex_solver::Bool
  robust::Bool
  temp_dir::String
  debug::Bool
  function Kraken(env, nmodes, mesh_density, clow, chigh, rmax, complex_solver, robust, temp_dir, debug)
    _check_env(Kraken, env)
    nmodes ≥ 1 || error("number of modes should be positive")
    mesh_density ≥ 0 || error("mesh density should be non-negative")
    clow ≥ 0.0 || error("clow should be non-negative")
    chigh > clow || error("chigh should be more than clow")
    rmax ≥ 0.0 || error("rmax should be non-negative")
    new{typeof(env)}(env, nmodes, mesh_density, clow, chigh, rmax, complex_solver, robust, temp_dir, debug)
  end
end

"""
    Kraken(env; kwargs...)

Create a Kraken propagation model.

Supported keyword arguments:
- `nmodes`: number of modes to use (default: 9999)
- `mesh_density`: number of mesh points per wavelength (default: 0, 0=auto)
- `clow`: lower limit of phase speed (default: 1300, 0=auto)
- `chigh`: upper limit of phase speed (default: 2500)
- `rmax`: largest range (in m) for field calculations (default: 0, 0=auto)
- `complex_solver`: use KrakenC for finding modes (default: true)
- `robust`: use robust (but slow) root finder (default: false)
- `temp_dir`: directory for temporary files (default: system temp directory)
- `debug`: debug mode (default: false)

Enabling debug mode will create a temporary directory with the Kraken input and output files.
This allows manual inspection of the files.
"""
function Kraken(env; nmodes=9999, mesh_density=0, clow=1300.0, chigh=2500.0, rmax=0.0, complex_solver=true, robust=false, temp_dir=tempdir(), debug=false)
  Kraken(env, nmodes, mesh_density, clow, chigh, rmax, complex_solver, robust, temp_dir, debug)
end

Base.show(io::IO, pm::Kraken) = print(io, "Kraken(⋯)")

### interface functions

function UnderwaterAcoustics.arrivals(pm::Kraken, tx1::AbstractAcousticSource, rx1::AbstractAcousticReceiver)
  mktempdir(pm.temp_dir; prefix="kraken_") do dirname
    # replace a single receiver with a grid of receivers at λ/10 spacing to sample the modes
    max_c = maximum(pm.env.soundspeed)
    λ = max_c / tx1.frequency
    D = maximum(pm.env.bathymetry)
    if pm.env.seabed isa MultilayerElasticBoundary
      # increase depth to include sediment layers
      D += sum(l -> l.h, pm.env.seabed.layers[1:end-1])
      max_c = 2 * max(max_c, maximum(l -> max(maximum(l.cₚ), maximum(l.cₛ)), pm.env.seabed.layers[1:end-1]))
    elseif pm.env.seabed isa ElasticBoundary
      max_c = 2 * max(max_c, maximum(pm.env.seabed.cₚ), maximum(pm.env.seabed.cₛ))
    else
      max_c = 2 * max(max_c, maximum(pm.env.seabed.c))
    end
    rxs = AcousticReceiverGrid2D(rx1.pos.x, range(-D, 0; length=ceil(Int, 10D/λ + 1)))
    _write_env(pm, tx1, rxs, dirname)
    _kraken(dirname, pm.complex_solver, pm.debug)
    ϕ, kᵣ, depths = _read_mod(pm, dirname)    # read mode shapes
    m, k, v = _read_grp(pm, dirname)          # read group velocity
    if all(0 .≤ v .≤ max_c)
      vs = SampledField(v; x=k)               # interpolate group velocity
      v = map(k -> vs(real(k)), kᵣ)
    else
      # replace invalid velocities by missing since KRAKEN gives zeros
      # and KRAKENC gives invalid values for multilayered elastic seabeds
      v = fill(missing, length(kᵣ))
    end
    ω = 2π * tx1.frequency
    map(1:min(length(kᵣ), pm.nmodes)) do i
      ψ = SampledField(ϕ[:,i]; z=-depths)
      ModeArrival{ComplexF64,typeof(ψ),Union{Missing,Float64},Float64}(i, kᵣ[i], ψ, v[i], ω/real(kᵣ[i]))
    end
  end
end

function UnderwaterAcoustics.acoustic_field(pm::Kraken, tx1::AbstractAcousticSource, rx::AcousticReceiverGrid2D; mode=:coherent)
  if mode === :coherent
    taskcode = 'C'
  elseif mode === :incoherent
    taskcode = 'I'
  else
    error("Unknown mode :" * string(mode))
  end
  fld = mktempdir(pm.temp_dir; prefix="kraken_") do dirname
    xrev, zrev = _write_env(pm, tx1, rx, dirname)
    _kraken(dirname, pm.complex_solver, pm.debug)
    _write_flp(pm, tx1, rx, dirname, mode)
    _field(dirname, pm.debug)
    _read_shd(joinpath(dirname, "model.shd"); xrev, zrev)
  end
  # conjugation required because Kraken uses the convention of cis(-kᵣR)
  conj.(fld .* db2amp(spl(tx1)))
end

function UnderwaterAcoustics.acoustic_field(pm::Kraken, tx1::AbstractAcousticSource, rx1::AbstractAcousticReceiver; mode=:coherent)
  if mode === :coherent
    taskcode = 'C'
  elseif mode === :incoherent
    taskcode = 'I'
  else
    error("Unknown mode :" * string(mode))
  end
  fld = mktempdir(pm.temp_dir; prefix="bellhop_") do dirname
    _write_env(pm, tx1, [rx1], dirname)
    _kraken(dirname, pm.complex_solver, pm.debug)
    _write_flp(pm, tx1, [rx1], dirname, mode)
    _field(dirname, pm.debug)
    _read_shd(joinpath(dirname, "model.shd"))[1]
  end
  # conjugation required because Kraken uses the convention of cis(-kᵣR)
  conj(fld * db2amp(spl(tx1)))
end

### helper functions

function _check_env(::Type{Kraken}, env)
  env.seabed isa FluidBoundary || env.seabed isa ElasticBoundary || env.seabed isa MultilayerElasticBoundary || error("seabed must be a FluidBoundary, ElasticBoundary or MultilayerElasticBoundary")
  env.surface isa FluidBoundary || error("surface must be a FluidBoundary")
  is_range_dependent(env.soundspeed) && error("Range-dependent soundspeed not supported")
  is_range_dependent(env.altimetry) && error("Range-dependent altimetry not supported")
  is_range_dependent(env.bathymetry) && error("Range-dependent bathymetry not supported")
  nothing
end

function _kraken(dirname, complex_solver, debug)
  outfilename = joinpath(dirname, "output.txt")
  try
    solver = complex_solver ? AcousticsToolbox_jll.krakenc() : AcousticsToolbox_jll.kraken()
    run(pipeline(ignorestatus(Cmd(`$solver model`; dir=dirname)); stdout=outfilename, stderr=outfilename))
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
  outfilename = joinpath(dirname, "output.txt")
  try
    run(pipeline(ignorestatus(Cmd(`$(AcousticsToolbox_jll.field()) model`; dir=dirname)); stdout=outfilename, stderr=outfilename))
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

function _write_flp(pm::Kraken, tx1, rx, dirname, mode)
  filename = joinpath(dirname, "model.flp")
  open(filename, "w") do io
    println(io, "/")                # take title from modes file
    print(io, "'RA ")               # point source, adiabatic modes
    print(io, mode == :incoherent ? 'I' : 'C')
    println(io, "'")
    @printf(io, "%i\n", pm.nmodes)  # maximum number of modes to use
    println(io, "1")                # number of profiles (for range dependence)
    println(io, "0.0")              # range (km) of first profile
    if length(rx) == 1
      _print_array(io, [location(rx[1]).x / 1000] )       # receiver ranges (km)
      _print_array(io, [-location(tx1).z])                # source depths (m)
      _print_array(io, [-location(rx[1]).z])              # receiver depths (m)
      _print_array(io, zeros(1))                          # receiver range displacements
    elseif rx isa AcousticReceiverGrid2D
      d = reverse(-rx.zrange)
      r = rx.xrange ./ 1000.0
      first(d) > last(d) && (d = reverse(d))
      first(r) > last(r) && (r = reverse(r))
      _print_array(io, r)                                 # receiver ranges (km)
      _print_array(io, [-location(tx1).z])                # source depths (m)
      _print_array(io, d)                                 # receiver depths (m)
      _print_array(io, zeros(length(d)))                  # receiver range displacements
    else
      error("Receivers must be on a 2D grid")
    end
  end
end

function _read_mod(pm::Kraken, dirname)
  filename = joinpath(dirname, "model.mod")
  open(filename, "r") do f
    reclen = 4 * read(f, UInt32)
    title = String(read!(f, zeros(UInt8, 80))) |> strip
    nfreq = read(f, UInt32) |> Int
    @assert nfreq == 1              # supports only one frequency
    nmedia = read(f, UInt32) |> Int
    ntot = read(f, UInt32) |> Int
    nmat = read(f, UInt32) |> Int
    seek(f, 4reclen)
    depths = read!(f, zeros(Float32, ntot))
    seek(f, 5reclen)
    nmodes = read(f, UInt32) |> Int
    rec = 6
    ϕ = zeros(ComplexF32, nmat, nmodes)
    for i in 1:nmodes
      seek(f, (rec + i) * reclen)
      x = read!(f, zeros(Float32, 2, nmat))
      ϕ[:,i] = complex.(x[1,:], x[2,:])
    end
    seek(f, (rec + nmodes + 1) * reclen)
    x = read!(f, zeros(Float32, 2, nmodes))
    kᵣ = complex.(x[1,:], x[2,:])
    (; ϕ, kᵣ, depths)
  end
end

function _read_grp(pm::Kraken, dirname)
  filename = joinpath(dirname, "model.prt")
  m = Int[]
  kᵣ = Float32[]
  v = Float32[]
  s = readlines(filename)
  i = findfirst(contains("Group Speed"), s)
  if i !== nothing
    i += 2
    while i <= length(s)
      f = split(strip(s[i]), r" +")
      length(f) == 5 || break
      push!(m, parse(Int, f[1]))
      push!(kᵣ, parse(Float32, f[2]))
      push!(v, parse(Float32, f[5]))
      i += 1
    end
  end
  (; m, kᵣ, v)
end
