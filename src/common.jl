struct ExecError <: Exception
  name::String
  details::Vector{String}
end

function Base.show(io::IO, e::ExecError)
  if length(e.details) == 1
    println(io, e.details[1])
  else
    println(io, "$(e.name) said:")
    for s ∈ e.details
      println(io, "  ", s)
    end
  end
end

function _check_err!(err, filename)
  output = false
  open(filename) do f
    for s in eachline(f)
      occursin("WARNING", uppercase(s)) && @warn s
      if output || occursin("ERROR", uppercase(s))
        push!(err, s)
        output = true
      end
    end
  end
end

function _write_env(pm, tx, rx, dirname; nbeams=0, taskcode=' ')
  all(location(tx1).x == 0 for tx1 ∈ tx) || error("Transmitters must be at (0, 0, z)")
  all(location(tx1).y == 0 for tx1 ∈ tx) || error("2D model requires transmitters in the x-z plane")
  all(location(rx1).x >= 0 for rx1 ∈ rx) || error("Receivers must be in the +x halfspace")
  all(location(rx1).y == 0 for rx1 ∈ rx) || error("2D model requires receivers in the x-z plane")
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
    nmedia = env.seabed isa MultilayerElasticBoundary ? length(env.seabed.layers) : 1
    println(io, nmedia)
    if length(rx) == 1
      maxr = location(only(rx)).x
    elseif rx isa AcousticReceiverGrid2D
      maxr = maximum(rx.xrange)
    else
      error("Receivers must be on a 2D grid")
    end
    ssp = env.soundspeed
    sspi = 'C'
    ssp isa SampledFieldZ && ssp.interp === :cubic && (sspi = 'S')
    surf = env.surface === RigidBoundary ? 'R' : env.surface === PressureReleaseBoundary ? 'V' : 'A'
    print(io, "'", sspi, surf, "WF")  # bottom attenuation in dB/wavelength, Francois-Garrison volume attenuation
    pm isa Kraken && pm.robust && print(io, ".")
    if pm isa Bellhop && is_range_dependent(env.altimetry)
      print(io, "*")
      _create_alt_bathy_file(joinpath(dirname, "model.ati"), env.altimetry, (q, p) -> -value(q, p), maxr, f)
    end
    println(io, "'")
    bathy = env.bathymetry
    waterdepth = maximum(bathy)
    @printf(io, "%0.1f %0.1f %0.1f %0.1f\n", env.temperature, env.salinity, env.pH, waterdepth/2)
    if surf == 'A'
      @printf(io, "0.0 %0.6f 0.0 %0.6f %0.6f /\n", env.surface.c, env.surface.ρ / env.density, in_dBperλ(env.surface.δ))
    end
    if pm isa Kraken
      λ = maximum(ssp) / f
      nmesh = ceil(Int, pm.nmesh_per_λ * λ)
      @printf(io, "%i %0.6f %0.6f\n", nmesh, env.surface.σ, waterdepth)
    else
      @printf(io, "0 0.0 %0.6f\n", waterdepth)
    end
    if is_constant(ssp)
      @printf(io, "0.0 %0.6f /\n", value(ssp))
      @printf(io, "%0.6f %0.6f /\n", waterdepth, value(ssp))
    elseif ssp isa SampledFieldZ
      zrange = sort!(vcat(collect(ssp.zrange), -waterdepth); rev=true)
      for z ∈ zrange
        @printf(io, "%0.6f %0.6f /\n", -z, ssp(z))
        z == -waterdepth && break
      end
    else
      for d ∈ range(0.0, waterdepth; length=_recommend_len(waterdepth, f))
        @printf(io, "%0.6f %0.6f /\n", d, ssp(-d))
      end
      floor(waterdepth) != waterdepth && @printf(io, "%0.6f %0.6f /\n", waterdepth, ssp(-waterdepth))
    end
    if nmedia > 1
      for l ∈ env.seabed.layers[1:end-1]
        ρ₁, ρ₂ = first(l.ρ), last(l.ρ)
        cₚ₁, cₚ₂ = first(l.cₚ), last(l.cₚ)
        cₛ₁, cₛ₂ = first(l.cₛ), last(l.cₛ)
        λ = max(cₚ₁, cₛ₁, cₚ₂, cₛ₂) / f
        nmesh = ceil(Int, 2 * pm.nmesh_per_λ * λ)    # Kraken manual recommends double the number of mesh points for elastic media
        @printf(io, "%i %0.6f %0.6f\n", nmesh, l.σ, waterdepth + l.h)
        @printf(io, "%0.6f %0.6f %0.6f %0.6f %0.6f %0.6f\n", waterdepth, cₚ₁, cₛ₁, ρ₁ / env.density, in_dBperλ(l.δₚ), in_dBperλ(l.δₛ))
        @printf(io, "%0.6f %0.6f %0.6f %0.6f /\n", waterdepth + l.h, cₚ₂, cₛ₂, ρ₂ / env.density)
        waterdepth += l.h
      end
    end
    seabed = env.seabed
    if seabed isa MultilayerElasticBoundary
      l = env.seabed.layers[end]
      seabed = ElasticBoundary(l.ρ, l.cₚ, l.cₛ, l.δₚ, l.δₛ)
    end
    if seabed isa ElasticBoundary && seabed.cₛ == 0
      seabed = FluidBoundary(seabed.ρ, seabed.cₚ, seabed.δₚ)
    end
    print(io, seabed === RigidBoundary ? "'R" : seabed === PressureReleaseBoundary ? "'V" : "'A")
    if is_range_dependent(bathy)
      print(io, "*")
      _create_alt_bathy_file(joinpath(dirname, "model.bty"), bathy, value, maxr, f)
    end
    @printf(io, "' %0.6f\n", seabed.σ)
    if seabed !== RigidBoundary && seabed !== PressureReleaseBoundary
      if seabed isa ElasticBoundary
        @printf(io, "%0.6f %0.6f %0.6f %0.6f %0.6f %0.6f\n", waterdepth, seabed.cₚ, seabed.cₛ, seabed.ρ / env.density, in_dBperλ(seabed.δₚ), in_dBperλ(seabed.δₛ))
      else
        @printf(io, "%0.6f %0.6f 0.0 %0.6f %0.6f /\n", waterdepth, seabed.c, seabed.ρ / env.density, in_dBperλ(seabed.δ))
      end
    end
    if pm isa Kraken
      @printf(io, "%0.6f  %0.6f\n", pm.clow, pm.chigh)
      @printf(io, "%0.6f\n", 1.01 * maxr / 1000.0)
    end
    _print_array(io, [-location(tx1).z for tx1 ∈ tx])
    if length(rx) == 1
      _print_array(io, [-location(rx[1]).z])
      pm isa Kraken || _print_array(io, [maxr / 1000.0])
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
      _print_array(io, d)
      pm isa Kraken || _print_array(io, r)
    else
      error("Receivers must be on a 2D grid")
    end
    if pm isa Bellhop
      bcode = pm.beam_type == :cartesian ? 'C' : pm.beam_type == :ray_centered ? 'R' : pm.beam_type == :gaussian ? 'B' : 'G'
      println(io, "'", taskcode, bcode, pm.beam_shift ? "S'" : "'")
      @printf(io, "%d\n", max(nbeams, pm.nbeams))
      @printf(io, "%0.6f %0.6f /\n", -rad2deg(pm.max_angle), -rad2deg(pm.min_angle))
      @printf(io, "0.0 %0.6f %0.6f\n", 1.01*waterdepth, 1.01 * maxr / 1000.0)
    end
  end
  xrev, zrev
end

function _print_array(io, a::AbstractVector)
  println(io, length(a))
  for a1 ∈ a
    @printf(io, "%0.6f ", a1)
  end
  println(io, "/")
end

function _recommend_len(x, f)
  # recommendation based on nominal half-wavelength spacing
  λ = 1500.0 / f
  clamp(round(Int, 2x / λ) + 1, 25, 1000)
end

function _create_alt_bathy_file(filename, data, func, maxr, f)
  open(filename, "w") do io
    interp = "L"
    if data isa SampledFieldX
      x = data.xrange
      data.interp !== :linear && (interp = "C")
    else
      x = range(0.0, maxr; length=_recommend_len(maxr, f))
    end
    println(io, "'", interp, "'")
    println(io, length(x))
    for i ∈ 1:length(x)
      @printf(io, "%0.6f %0.6f\n", x[i]/1000.0, func(data, (x[i], 0.0)))
    end
  end
end

function _read_rays(filename)
  rays = RayArrival{Float64,Float64}[]
  open(filename, "r") do io
    [readline(io) for i ∈ 1:7]
    while !eof(io)
      s = strip(readline(io))
      length(s) == 0 && break
      aod = parse(Float64, s)
      pts, sb, bb = parse.(Int, split(strip(readline(io)), r" +"))
      raypath = Vector{XYZ{NTuple{3,Float64}}}(undef, pts)
      for k ∈ 1:pts
        eof(io) && break
        x, d = parse.(Float64, split(strip(readline(io)), r" +"))
        raypath[k] = (x=x, y=0.0, z=-d)
      end
      push!(rays, RayArrival(NaN64, complex(NaN64, NaN64), sb, bb, -deg2rad(aod), NaN64, raypath))
    end
  end
  rays
end

function _read_arr(filename)
  T = XYZ{NTuple{3,Float64}}
  arrivals = RayArrival{Float64,Float64,Float64,Float64,Vector{T}}[]
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
            (A == 0 || isnan(A)) && continue
            push!(arrivals, eltype(arrivals)(tr, A * cis(deg2rad(ph) - 2π * f * complex(0, ti)), sb, bb, -deg2rad(aod), deg2rad(aoa), T[]))
          end
        end
      end
    end
  end
  sort(arrivals; by = a -> a.t)
end

function _read_shd(filename; xrev=false, zrev=false)
  open(filename, "r") do io
    r = read(io, UInt32)
    seek(io, 4r)
    b = Array{UInt8}(undef, 10)
    read!(io, b)  # strip(String(b)) == "rectilin" for Bellhop, but blank for Kraken
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
