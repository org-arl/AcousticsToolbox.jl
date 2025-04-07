export Bellhop

const BELLHOP = Ref{Cmd}(AcousticsToolbox_jll.bellhop())

"""
A propagation model based on the FORTRAN OALIB Bellhop model.
"""
struct Bellhop{T} <: AbstractRayPropagationModel
  env::T
  nbeams::Int
  minangle::Float32
  maxangle::Float32
  gaussian::Bool
  debug::Bool
  function Bellhop(env, nbeams, minangle, maxangle, gaussian, debug)
    checkenv(env)
    nbeams < 0 && (nbeams = 0)
    -π/2 ≤ minangle ≤ π/2 || throw(ArgumentError("minangle should be between -π/2 and π/2"))
    -π/2 ≤ maxangle ≤ π/2 || throw(ArgumentError("maxangle should be between -π/2 and π/2"))
    minangle < maxangle || throw(ArgumentError("maxangle should be more than minangle"))
    new{typeof(env)}(env, nbeams, Float32(minangle), Float32(maxangle), gaussian, debug)
  end
end

"""
    Bellhop(env; gaussian=false, debug=false)
    Bellhop(env, nbeams, minangle, maxangle, gaussian, debug)

Create a Bellhop propagation model.
"""
Bellhop(env; gaussian=false, debug=false) = Bellhop(env, 0, -deg2rad(80), deg2rad(80), gaussian, debug)

### interface functions

function UnderwaterAcoustics.arrivals(pm::Bellhop, tx1::AbstractAcousticSource, rx1::AbstractAcousticReceiver; paths=true)
  mktempdir(prefix="bellhop_") do dirname
    writeenv(pm, [tx1], [rx1], "A", dirname)
    bellhop(dirname, pm.debug)
    arr = readarrivals(joinpath(dirname, "model.arr"))
    if paths
      writeenv(pm, [tx1], [rx1], "E", dirname)
      bellhop(dirname, pm.debug)
      arr2 = readrays(joinpath(dirname, "model.ray"))
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
    taskcode = "C"
  elseif mode === :incoherent
    taskcode = "I"
  elseif mode === :semicoherent
    taskcode = "S"
  else
    throw(ArgumentError("Unknown mode :" * string(mode)))
  end
  mktempdir(prefix="bellhop_") do dirname
    xrev, zrev = writeenv(pm, [tx1], rx, taskcode, dirname)
    bellhop(dirname, pm.debug)
    readshd(joinpath(dirname, "model.shd"); xrev=xrev, zrev=zrev)
  end
end

function UnderwaterAcoustics.acoustic_field(pm::Bellhop, tx1::AbstractAcousticSource, rx1::AbstractAcousticReceiver; mode=:coherent)
  if mode === :coherent
    taskcode = "C"
  elseif mode === :incoherent
    taskcode = "I"
  elseif mode === :semicoherent
    taskcode = "S"
  else
    throw(ArgumentError("Unknown mode :" * string(mode)))
  end
  mktempdir(prefix="bellhop_") do dirname
    writeenv(pm, [tx1], [rx1], taskcode, dirname)
    bellhop(dirname, pm.debug)
    readshd(joinpath(dirname, "model.shd"))[1]
  end
end

### helper functions

struct BellhopError <: Exception
  details::Vector{String}
end

function Base.show(io::IO, e::BellhopError)
  if length(e.details) == 1
    println(io, e.details[1])
  else
    println(io, "Bellhop said:")
    for s ∈ e.details
      println(io, "  ", s)
    end
  end
end

function bellhop(dirname, debug)
  infilebase = joinpath(dirname, "model")
  outfilename = joinpath(dirname, "output.txt")
  try
    run(pipeline(ignorestatus(`$(BELLHOP[]) $infilebase`); stdout=outfilename, stderr=outfilename))
    if debug
      @info "Bellhop run completed in $dirname, press ENTER to delete intermediate files..."
      readline()
    end
  catch
    throw(BellhopError(["Unable to execute bellhop.exe"]))
  end
  err = String[]
  checkerr!(err, outfilename)
  checkerr!(err, joinpath(dirname, "model.prt"))
  if length(err) > 0
    throw(BellhopError(err))
  end
end

function checkerr!(err, filename)
  output = false
  open(filename) do f
    for s in eachline(f)
      if output || occursin("ERROR", uppercase(s))
        push!(err, s)
        output = true
      end
    end
  end
end

function checkenv(env)
  env.seabed isa FluidBoundary || throw(ArgumentError("seabed must be a FluidBoundary"))
  env.surface === PressureReleaseBoundary || throw(ArgumentError("surface must be a PressureReleaseBoundary"))
  is_range_dependent(env.soundspeed) && throw(ArgumentError("range-dependent soundspeed not supported"))
  mktempdir(prefix="bellhop_") do dirname
    try
      bellhop(dirname, false)
    catch e
      e isa BellhopError && e.details == ["Unable to execute bellhop.exe"] && throw(e)
    end
  end
  nothing
end

function writeenv(pm::Bellhop, tx, rx, taskcode, dirname; minangle=pm.minangle, maxangle=pm.maxangle, nbeams=pm.nbeams)
  all(location(tx1).x == 0 for tx1 ∈ tx) || throw(ArgumentError("Bellhop requires transmitters at (0, 0, z)"))
  all(location(tx1).y == 0 for tx1 ∈ tx) || throw(ArgumentError("Bellhop 2D requires transmitters in the x-z plane"))
  all(location(rx1).x >= 0 for rx1 ∈ rx) || throw(ArgumentError("Bellhop requires receivers to be in the +x halfspace"))
  all(location(rx1).y == 0 for rx1 ∈ rx) || throw(ArgumentError("Bellhop 2D requires receivers in the x-z plane"))
  xrev = false
  zrev = false
  env = pm.env
  name = split(basename(dirname), "_")[end]
  filename = joinpath(dirname, "model.env")
  open(filename, "w") do io
    println(io, "'", name, "'")
    flist = [frequency(tx1) for tx1 ∈ tx]
    f = sum(flist) / length(flist)
    maximum(abs, flist .- f) / f > 0.2 && @warn("Source frequency varies by more than 20% from nominal frequency")
    @printf(io, "%0.6f\n", f)
    println(io, "1")
    if length(rx) == 1
      maxr = location(only(rx)).x
    elseif rx isa AcousticReceiverGrid2D
      maxr = maximum(rx.xrange)
    else
      throw(ArgumentError("Receivers must be on a 2D grid"))
    end
    ssp = env.soundspeed
    sspi = "S"
    ssp isa SampledFieldZ && ssp.interp === :linear && (sspi = "C")
    print(io, "'", sspi, "VWT")
    if is_range_dependent(env.altimetry)
      print(io, "*")
      createadfile(joinpath(dirname, "model.ati"), env.altimetry, (q, p) -> -value(q, p), maxr, f)
    end
    println(io, "'")
    bathy = env.bathymetry
    if is_constant(bathy)
      waterdepth = value(bathy)
    else
      waterdepth = maximum(x -> value(bathy, (x, 0, 0)), range(0.0, maxr; length=recommendlength(maxr, f)))
    end
    @printf(io, "1 0.0 %0.6f\n", waterdepth)
    if is_constant(ssp)
      @printf(io, "0.0 %0.6f /\n", value(ssp))
      @printf(io, "%0.6f %0.6f /\n", waterdepth, value(ssp))
    elseif ssp isa SampledFieldZ
      for z ∈ ssp.zrange
        @printf(io, "%0.6f %0.6f /\n", -z, ssp(z))
      end
    else
      for d ∈ range(0.0, waterdepth; length=recommendlength(waterdepth, f))
        @printf(io, "%0.6f %0.6f /\n", d, ssp(-d))
      end
      floor(waterdepth) != waterdepth && @printf(io, "%0.6f %0.6f /\n", waterdepth, ssp(-waterdepth))
    end
    print(io, "'A")
    if is_range_dependent(bathy)
      print(io, "*")
      createadfile(joinpath(dirname, "model.bty"), bathy, value, maxr, f)
    end
    println(io, "' 0.0") # bottom roughness = 0
    @printf(io, "%0.6f %0.6f 0.0 %0.6f %0.6f /\n", waterdepth, env.seabed.c, env.seabed.ρ / 1000, in_dBperλ(env.seabed.δ))
    printarray(io, [-location(tx1).z for tx1 ∈ tx])
    if length(rx) == 1
      printarray(io, [-location(rx[1]).z])
      printarray(io, [maxr / 1000.0])
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
      printarray(io, d)
      printarray(io, r)
    end
    println(io, "'", taskcode, pm.gaussian ? "B'" : "'")
    @printf(io, "%d\n", nbeams)
    @printf(io, "%0.6f %0.6f /\n", rad2deg(minangle), rad2deg(maxangle))
    @printf(io, "0.0 %0.6f %0.6f\n", 1.01*waterdepth, 1.01 * maxr / 1000.0)
  end
  xrev, zrev
end

function printarray(io, a::AbstractVector)
  println(io, length(a))
  for a1 ∈ a
    @printf(io, "%0.6f ", a1)
  end
  println(io, "/")
end

function recommendlength(x, f)
  # recommendation based on nominal half-wavelength spacing
  λ = 1500.0 / f
  clamp(round(Int, 2x / λ) + 1, 25, 1000)
end

function createadfile(filename, data, func, maxr, f)
  open(filename, "w") do io
    interp = "L"
    if data isa SampledDepth || data isa SampledAltitude
      x = data.x
      data.interp !== :linear && (interp = "C")
    else
      x = range(0.0, maxr; length=recommendlength(maxr, f))
    end
    println(io, "'", interp, "'")
    println(io, length(x))
    for i ∈ 1:length(x)
      @printf(io, "%0.6f %0.6f\n", x[i]/1000.0, func(data, x[i], 0.0))
    end
  end
end

function readrays(filename)
  rays = RayArrival{Float64,Float64}[]
  open(filename, "r") do io
    [readline(io) for i ∈ 1:7]
    while !eof(io)
      s = strip(readline(io))
      length(s) == 0 && break
      aod = parse(Float64, s)
      pts, sb, bb = parse.(Int, split(strip(readline(io)), r" +"))
      raypath = Array{NTuple{3,Float64}}(undef, pts)
      for k ∈ 1:pts
        eof(io) && break
        x, d = parse.(Float64, split(strip(readline(io)), r" +"))
        raypath[k] = (x, 0.0, -d)
      end
      push!(rays, RayArrival(NaN64, complex(NaN64, NaN64), sb, bb, -deg2rad(aod), NaN64, raypath))
    end
  end
  rays
end

function readarrivals(filename)
  arrivals = RayArrival{Float64,Float64,Float64,Float64,Union{Missing,Vector{NTuple{3,Float64}}}}[]
  open(filename, "r") do io
    s = strip(readline(io))
    if occursin("2D", s)
      f = parse(Float64, strip(readline(io)))
      v = split(strip(readline(io)) ,r" +")
      n = parse(Int, v[1])
      txdepth = parse.(Float64, v[2:end])
      n == length(txdepth) || error("Wrong number of txdepth entries in arrivals")
      v = split(strip(readline(io)) ,r" +")
      n = parse(Int, v[1])
      rxdepth = parse.(Float64, v[2:end])
      n == length(rxdepth) || error("Wrong number of rxdepth entries in arrivals")
      v = split(strip(readline(io)) ,r" +")
      n = parse(Int, v[1])
      rxrange = parse.(Float64, v[2:end])
      n == length(rxrange) || error("Wrong number of rxrange entries in arrivals")
    else
      v = split(s ,r" +")
      f = parse(Float64, v[1])
      n1, n2, n3 = parse.(Int, v[2:4])
      txdepth = parse.(Float64, split(strip(readline(io)) ,r" +"))
      rxdepth = parse.(Float64, split(strip(readline(io)) ,r" +"))
      rxrange = parse.(Float64, split(strip(readline(io)) ,r" +"))
      n1 == length(txdepth) || error("Wrong number of txdepth entries in arrivals")
      n2 == length(rxdepth) || error("Wrong number of rxdepth entries in arrivals")
      n3 == length(rxrange) || error("Wrong number of rxrange entries in arrivals")
    end
    for j ∈ 1:length(txdepth)
      readline(io)
      for k ∈ 1:length(rxdepth)
        for m ∈ 1:length(rxrange)
          count = parse(Int, strip(readline(io)))
          for n ∈ 1:count
            v = split(strip(readline(io)) ,r" +")
            length(v) == 8 || error("Wrong number of data entries in arrivals")
            A, ph, tr, ti, aod, aoa = parse.(Float64, v[1:6])
            sb, bb = parse.(Int, v[7:8])
            push!(arrivals, eltype(arrivals)(tr, A * cis(deg2rad(ph) - 2π * f * complex(0, ti)), sb, bb, -deg2rad(aod), deg2rad(aoa), missing))
          end
        end
      end
    end
  end
  sort(arrivals; by = a -> a.t)
end

function readshd(filename; xrev=false, zrev=false)
  open(filename, "r") do io
    r = read(io, UInt32)
    seek(io, 4r)
    b = Array{UInt8}(undef, 10)
    read!(io, b)
    strip(String(b)) == "rectilin" || error("Bad shd file format: incorrect ptype")
    seek(io, 8r)
    nfreq = read(io, UInt32)
    nfreq == 1 || error("Bad shd file format: incorrect nfreq")
    nθ = read(io, UInt32)
    nθ == 1 || error("Bad shd file format: incorrect nθ")
    nsx = read(io, UInt32)
    nsy = read(io, UInt32)
    nsd = read(io, UInt32)
    nsd == 1 || error("Bad shd file format: incorrect nsd")
    nrd = read(io, UInt32)
    nrr = read(io, UInt32)
    pressure = Array{ComplexF32}(undef, nrr, nrd)
    for ird ∈ 0:nrd-1
      recnum = 10 + ird
      seek(io, recnum * 4r)
      temp = Array{ComplexF32}(undef, nrr)
      read!(io, temp)
      pressure[:,ird+1] .= -temp    # negative because Bellhop seems to have a 180° phase inversion
    end
    xrev && (pressure = reverse(pressure; dims=1))
    zrev || (pressure = reverse(pressure; dims=2))
    pressure
  end
end
