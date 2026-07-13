using TestItems

# Tests for the KrakenJL solver (src/krakenjl/, GPL-3.0-or-later).
# Golden reference constants below were generated from the Fortran Kraken
# (AcousticsToolbox_jll, OALIB 2024_12_25 build) via the Kraken wrapper in
# this package.

@testitem "krakenjl-models" begin
  using UnderwaterAcoustics
  m = models()
  @test KrakenJL ∈ m
  env = UnderwaterEnvironment(bathymetry = 20.0u"m")
  pm = KrakenJL(env)
  @test pm isa KrakenJL
  @test_throws ErrorException KrakenJL(env; nmodes = 0)
  @test_throws ErrorException KrakenJL(env; mesh_density = -1)
  @test_throws ErrorException KrakenJL(env; clow = 2000.0, chigh = 1500.0)
  @test_throws ErrorException KrakenJL(env; threads = 0)
end

@testitem "krakenjl-pekeris-modes" begin
  using UnderwaterAcoustics
  # mirrors test_kraken.jl "kraken-pekeris-modes"
  env = UnderwaterEnvironment(
    bathymetry = 5000.0,
    soundspeed = 1500.0,
    density = 1000.0,
    seabed = FluidBoundary(2000.0, 2000.0)
  )
  tx = AcousticSource(0.0, -500.0, 10.0)
  rx = AcousticReceiver(200000.0, -2500.0)
  rxs = AcousticReceiverGrid2D(200000.0:10.0:220000.0, -2500.0)
  pm = PekerisModeSolver(env)
  arr = arrivals(pm, tx, rx)
  x = transmission_loss(pm, tx, rxs)
  for complex_solver ∈ [false, true]
    pm1 = KrakenJL(env; chigh=2000.0, complex_solver)
    arr1 = arrivals(pm1, tx, rx)
    @test arr1 isa AbstractArray{<:UnderwaterAcoustics.ModeArrival}
    @test length(arr1) == length(arr)
    @test [a.kᵣ for a ∈ arr1] ≈ [a.kᵣ for a ∈ arr] atol=0.001
    x1 = transmission_loss(pm1, tx, rxs)
    @test size(x1) == size(x)
    @test sum(abs, x1 .- x) / length(x) < 1
  end
  # KRAKENC group speeds should be physical; KRAKEN (real) reports zeros
  arrc = arrivals(KrakenJL(env; chigh=2000.0), tx, rx)
  @test all(0 .< [a.v for a ∈ arrc] .< 1500.0)
  arrr = arrivals(KrakenJL(env; chigh=2000.0, complex_solver=false), tx, rx)
  @test all(iszero, [a.v for a ∈ arrr])
end

@testitem "krakenjl-golden" begin
  using UnderwaterAcoustics
  # golden (Fortran Kraken, complex_solver=true) references
  env = UnderwaterEnvironment(bathymetry = 5000.0, soundspeed = 1500.0,
    density = 1000.0, seabed = FluidBoundary(2000.0, 2000.0))
  tx = AcousticSource(0.0, -500.0, 10.0)
  arr = arrivals(KrakenJL(env; chigh=2000.0), tx, AcousticReceiver(200000.0, -2500.0))
  @test length(arr) == 44
  @test real(arr[1].kᵣ) ≈ 0.041883323 atol=1e-7
  @test real(arr[2].kᵣ) ≈ 0.041869581 atol=1e-7
  @test real(arr[10].kᵣ) ≈ 0.041426945 atol=1e-7
  @test real(arr[44].kᵣ) ≈ 0.031728230 atol=1e-7
  @test imag(arr[1].kᵣ) ≈ -7.119e-10 rtol=0.01
  rxs = AcousticReceiverGrid2D(200000.0:10.0:220000.0, -2500.0)
  x = transmission_loss(KrakenJL(env; chigh=2000.0), tx, rxs)
  @test x[1] ≈ 93.2741 atol=0.1
  @test x[500] ≈ 86.6302 atol=0.1
  @test x[2001] ≈ 95.2300 atol=0.1
  # lossy elastic seabed (KRAKENC only; the real KRAKEN caps chigh at the
  # shear speed and finds no modes, as the Fortran does)
  env2 = UnderwaterEnvironment(bathymetry = 100.0, soundspeed = 1500.0,
    seabed = ElasticBoundary(2000.0, 1800.0, 400.0, 0.1, 0.2))
  tx2 = AcousticSource(0.0, -25.0, 300.0)
  arr2 = arrivals(KrakenJL(env2), tx2, AcousticReceiver(3000.0, -60.0))
  @test length(arr2) == 32
  @test real(arr2[1].kᵣ) ≈ 1.256259322 atol=1e-6
  @test imag(arr2[1].kᵣ) ≈ -6.218312e-6 rtol=0.01
  rxs2 = AcousticReceiverGrid2D(1000.0:10.0:3000.0, -60.0)
  x2 = transmission_loss(KrakenJL(env2), tx2, rxs2)
  @test x2[1] ≈ 61.1224 atol=0.1
  @test x2[100] ≈ 64.9941 atol=0.1
  @test x2[201] ≈ 58.1588 atol=0.1
end

@testitem "krakenjl-vs-kraken" begin
  using UnderwaterAcoustics
  using Statistics
  # multilayer elastic seabed and spline SSP, cross-checked at runtime
  env1 = UnderwaterEnvironment(bathymetry = 100.0, soundspeed = 1500.0,
    seabed = MultilayerElasticBoundary([
      (h=20.0, ρ=1500.0, cₚ=1600.0, cₛ=200.0, δₚ=0.1, δₛ=0.2),
      (h=Inf, ρ=2000.0, cₚ=1800.0, cₛ=400.0, δₚ=0.1, δₛ=0.2)]))
  tx1 = AcousticSource(0.0, -25.0, 300.0)
  rxs1 = AcousticReceiverGrid2D(1000.0:10.0:3000.0, -60.0)
  munk(z) = 1500 * (1 + 0.00737 * (2 * (z - 1300) / 1300 - 1 + exp(-2 * (z - 1300) / 1300)))
  zvals = 0.0:250.0:5000.0
  env2 = UnderwaterEnvironment(bathymetry = 5000.0,
    soundspeed = SampledField(munk.(zvals); z=-zvals, interp=CubicSpline()),
    seabed = FluidBoundary(1900.0, 1650.0, 0.1))
  tx2 = AcousticSource(0.0, -1000.0, 50.0)
  rxs2 = AcousticReceiverGrid2D(10000.0:100.0:60000.0, -1000.0)
  for (env, tx, rxs) ∈ [(env1, tx1, rxs1), (env2, tx2, rxs2)]
    x1 = transmission_loss(KrakenJL(env), tx, rxs)
    x2 = transmission_loss(Kraken(env), tx, rxs)
    @test median(abs.(x1 .- x2)) < 0.1
    y1 = transmission_loss(KrakenJL(env), tx, rxs; mode=:incoherent)
    y2 = transmission_loss(Kraken(env), tx, rxs; mode=:incoherent)
    @test median(abs.(y1 .- y2)) < 0.1
  end
end

@testitem "krakenjl-threads" begin
  using UnderwaterAcoustics
  # serial vs threaded results must be identical: each mode / range column is
  # written by exactly one task, so there is no reduction reordering
  env = UnderwaterEnvironment(bathymetry = 5000.0, soundspeed = 1500.0,
    density = 1000.0, seabed = FluidBoundary(2000.0, 2000.0))
  tx = AcousticSource(0.0, -500.0, 10.0)
  rxs = AcousticReceiverGrid2D(200000.0:100.0:220000.0, -4900.0:100.0:-100.0)
  x1 = acoustic_field(KrakenJL(env; chigh=2000.0, threads=1), tx, rxs)
  x4 = acoustic_field(KrakenJL(env; chigh=2000.0, threads=4), tx, rxs)
  @test x1 == x4
  @test x4 == acoustic_field(KrakenJL(env; chigh=2000.0, threads=4), tx, rxs)
  a1 = arrivals(KrakenJL(env; chigh=2000.0, threads=1), tx, AcousticReceiver(200000.0, -2500.0))
  a4 = arrivals(KrakenJL(env; chigh=2000.0, threads=4), tx, AcousticReceiver(200000.0, -2500.0))
  @test [a.kᵣ for a ∈ a1] == [a.kᵣ for a ∈ a4]
  @test all([a1[i].ψ(-100.0) == a4[i].ψ(-100.0) for i ∈ eachindex(a1)])
end
