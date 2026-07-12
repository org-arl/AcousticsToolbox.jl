using TestItems

# Tests for the BellhopJL solver (src/bellhopjl/, GPL-3.0-or-later).
# Golden reference constants below were generated from the Fortran Bellhop
# (AcousticsToolbox_jll, OALIB build) via the Bellhop wrapper in this package.

@testitem "bellhopjl-models" begin
  using UnderwaterAcoustics
  m = models()
  @test BellhopJL ∈ m
  env = UnderwaterEnvironment(bathymetry = 20.0u"m")
  pm = BellhopJL(env)
  @test pm isa BellhopJL
  @test_throws ErrorException BellhopJL(env; min_angle = 0.0, max_angle = 0.0)
  @test_throws ErrorException BellhopJL(env; beam_type = :cerveny)
end

@testitem "bellhopjl-arrivals" begin
  using UnderwaterAcoustics
  using UnderwaterAcoustics: distance
  # Pekeris golden scenario: 100 m deep, sand seabed, 500 Hz, tx at 25 m
  env = UnderwaterEnvironment(bathymetry = 100.0, soundspeed = 1500.0,
    seabed = FluidBoundary(1900.0, 1650.0, 0.0))
  pm = BellhopJL(env)
  tx = AcousticSource(0.0, -25.0, 500.0)
  rx = AcousticReceiver(1000.0, -60.0)
  arr = arrivals(pm, tx, rx)
  @test arr isa AbstractArray{<:UnderwaterAcoustics.RayArrival}
  @test length(arr) ≥ 10
  @test all([arr[j].t ≥ arr[j-1].t for j ∈ 2:length(arr)])
  # golden (Fortran) reference: first arrivals at r=1000, z=-60
  @test arr[1].t ≈ 0.66707480 atol=1e-5
  @test arr[2].t ≈ 0.66907066 atol=1e-5
  @test arr[3].t ≈ 0.67106044 atol=1e-5
  @test [(arr[j].ns, arr[j].nb) for j ∈ 1:4] == [(0,0), (1,0), (0,1), (1,1)]
  @test abs(arr[1].ϕ) ≈ 9.9763e-4 rtol=0.01     # direct path, ≈ 1/r
  @test real(arr[2].ϕ) < 0                      # surface bounce flips phase
  @test abs(arr[3].ϕ) ≈ hypot(-5.42729563e-4, 8.29995896e-4) rtol=0.02
  @test arr[1].θₛ ≈ -0.034986 atol=1e-4
  @test arr[1].θᵣ ≈ 0.034986 atol=1e-4
  # eigenray paths
  @test all([a.path[1] == (x=0.0, y=0.0, z=-25.0) for a ∈ arr])
  @test all([distance(a.path[end], (x=1000.0, y=0.0, z=-60.0)) < 20.0 for a ∈ arr[1:4]])
end

@testitem "bellhopjl-field" begin
  using UnderwaterAcoustics
  # deep isovelocity water, single direct path: |p| ≈ 1/r
  env = UnderwaterEnvironment(bathymetry = 5000.0, soundspeed = 1500.0,
    seabed = FluidBoundary(water_density(), 1500.0, 0.0))
  pm = BellhopJL(env; min_angle = -30°, max_angle = 30°)
  tx = AcousticSource((x=0.0, z=-2500.0), 300.0)
  x = acoustic_field(pm, tx, AcousticReceiver((x=1000.0, z=-2500.0)))
  @test x isa Complex
  @test abs(x) ≈ 1/1000 rtol=0.02
  x′ = acoustic_field(pm, tx, AcousticReceiver((x=1000.0, z=-2500.0)); mode=:incoherent)
  @test abs(x′) ≈ 1/1000 rtol=0.02
  y = transmission_loss(pm, tx, AcousticReceiver((x=1000.0, z=-2500.0)))
  @test -20 * log10(abs(x)) ≈ y atol=0.1
end

@testitem "bellhopjl-pekeris-tl" begin
  using UnderwaterAcoustics
  # golden (Fortran Bellhop) TL values on the Pekeris scenario grid
  env = UnderwaterEnvironment(bathymetry = 100.0, soundspeed = 1500.0,
    seabed = FluidBoundary(1900.0, 1650.0, 0.0))
  tx = AcousticSource(0.0, -25.0, 500.0)
  rxs = AcousticReceiverGrid2D(range(200.0, 5000.0; length=120),
                               range(-95.0, -5.0; length=45))
  pm = BellhopJL(env)
  xloss = transmission_loss(pm, tx, rxs)
  @test size(xloss) == (120, 45)
  @test xloss[10, 10] ≈ 49.1411 atol=0.1     # r=563 m,  z=-76.6 m
  @test xloss[50, 20] ≈ 53.5504 atol=0.1     # r=2176 m, z=-56.1 m
  @test xloss[100, 40] ≈ 56.6454 atol=0.1    # r=4193 m, z=-15.2 m
  xloss′ = transmission_loss(pm, tx, rxs; mode=:incoherent)
  @test xloss′[10, 10] ≈ 47.8685 atol=0.1
  @test xloss′[50, 20] ≈ 54.0370 atol=0.1
  @test xloss′[30, 5] ≈ 52.1272 atol=0.1
end

@testitem "bellhopjl-munk-tl" begin
  using UnderwaterAcoustics
  # Munk spline-SSP profile, Gaussian beams, incoherent TL vs Fortran golden
  munk(z) = 1500 * (1 + 0.00737 * (2 * (z - 1300) / 1300 - 1 + exp(-2 * (z - 1300) / 1300)))
  zvals = 0.0:250.0:5000.0
  ssp = SampledField(munk.(zvals); z=-zvals, interp=CubicSpline())
  env = UnderwaterEnvironment(bathymetry = 5000.0, soundspeed = ssp,
    seabed = FluidBoundary(1900.0, 1650.0, 0.0))
  tx = AcousticSource(0.0, -1000.0, 50.0)
  rxs = AcousticReceiverGrid2D(range(1000.0, 60_000.0; length=150),
                               range(-4900.0, -100.0; length=49))
  pm = BellhopJL(env; beam_type = :gaussian)
  xloss = transmission_loss(pm, tx, rxs; mode=:incoherent)
  @test size(xloss) == (150, 49)
  @test xloss[10, 10] ≈ 72.0804 atol=0.5     # r=4.6 km,  z=-4000 m
  @test xloss[50, 20] ≈ 80.3815 atol=0.5     # r=20.4 km, z=-3000 m
  @test xloss[100, 40] ≈ 80.4499 atol=0.5    # r=40.2 km, z=-1000 m
  @test xloss[30, 5] ≈ 77.0737 atol=0.5      # r=12.5 km, z=-4500 m
end

@testitem "bellhopjl-analytic-ray" begin
  # linear soundspeed gradient ⇒ circular-arc rays (internal API)
  using AcousticsToolbox.BellhopJLCore: CLinearSSP, Env2D, BeamParams,
    flat_boundary, VacuumBC, RigidBC, trace_ray
  c0, g, depth = 1500.0, 0.05, 5000.0
  ssp = CLinearSSP([0.0, depth], [c0, c0 + g * depth])
  env = Env2D(100.0, ssp, flat_boundary(0.0; top=true),
              flat_boundary(depth; top=false), VacuumBC(), RigidBC())
  beam = BeamParams(deltas=10.0, box_r=3000.0, box_z=depth, nbeams=1,
                    αmin=0.0, αmax=0.0)
  z0, θ0 = 1000.0, deg2rad(-5)
  ray = trace_ray(env, θ0, z0, beam)
  R = (c0 + g * z0) / (g * cos(θ0))
  xc, zc = R * sin(θ0), z0 - R * cos(θ0)
  for pt ∈ ray[2:10:end]
    @test hypot(pt.x[1] - xc, pt.x[2] - zc) ≈ R rtol=2e-4
  end
  # isovelocity: exact travel time
  env2 = Env2D(100.0, CLinearSSP([0.0, depth], [c0, c0]),
               flat_boundary(0.0; top=true), flat_boundary(depth; top=false),
               VacuumBC(), RigidBC())
  ray2 = trace_ray(env2, deg2rad(10), 2500.0, beam)
  for pt ∈ ray2
    s = hypot(pt.x[1], pt.x[2] - 2500.0)
    @test real(pt.τ) ≈ s / c0 atol=1e-12
  end
end

@testitem "∂bellhopjl" begin
  using UnderwaterAcoustics
  using DifferentiationInterface
  import ForwardDiff, FiniteDifferences
  fd = AutoFiniteDifferences(fdm=FiniteDifferences.central_fdm(5, 1))
  function ℳ((D, R, d1, d2, f, c))
    env = UnderwaterEnvironment(bathymetry = D, soundspeed = c, seabed = SandySilt)
    pm = BellhopJL(env; beam_type = :gaussian, nbeams = 100)
    transmission_loss(pm, AcousticSource((x=0.0, z=-d1), f), AcousticReceiver((x=R, z=-d2)))
  end
  x = [20.0, 100.0, 5.0, 10.0, 5000.0, 1500.0]
  @test gradient(ℳ, AutoForwardDiff(), x) ≈ gradient(ℳ, fd, x) atol=1e-3
end

@testitem "bellhopjl-vs-bellhop" begin
  using UnderwaterAcoustics
  using Statistics
  # both models live in this package: incoherent TL on a Pekeris grid must
  # agree closely (amplitude-only comparison, robust to stepping differences)
  env = UnderwaterEnvironment(bathymetry = 100.0, soundspeed = 1500.0,
    seabed = FluidBoundary(1900.0, 1650.0, 0.0))
  tx = AcousticSource(0.0, -25.0, 500.0)
  rxs = AcousticReceiverGrid2D(range(200.0, 5000.0; length=60),
                               range(-95.0, -5.0; length=23))
  tl1 = transmission_loss(BellhopJL(env), tx, rxs; mode=:incoherent)
  tl2 = transmission_loss(Bellhop(env), tx, rxs; mode=:incoherent)
  @test size(tl1) == size(tl2)
  @test median(abs.(tl1 .- tl2)) < 0.1
  # coherent fields should also agree tightly on this scenario
  tl1c = transmission_loss(BellhopJL(env), tx, rxs)
  tl2c = transmission_loss(Bellhop(env), tx, rxs)
  @test median(abs.(tl1c .- tl2c)) < 0.1
end

@testitem "bellhopjl-threads" begin
  using UnderwaterAcoustics
  # serial vs threaded results must agree (bitwise-reproducible per thread
  # count; serial vs threaded differ only by floating-point reassociation)
  env = UnderwaterEnvironment(bathymetry = 100.0, soundspeed = 1500.0,
    seabed = FluidBoundary(1900.0, 1650.0, 0.0))
  tx = AcousticSource(0.0, -25.0, 500.0)
  rxs = AcousticReceiverGrid2D(range(200.0, 5000.0; length=120),
                               range(-95.0, -5.0; length=45))
  pm1 = BellhopJL(env; threads=1)
  pm4 = BellhopJL(env; threads=4)
  @test_throws ErrorException BellhopJL(env; threads=0)
  x1 = acoustic_field(pm1, tx, rxs)
  x4 = acoustic_field(pm4, tx, rxs)
  @test x4 == acoustic_field(pm4, tx, rxs)       # deterministic
  @test isapprox(x1, x4; rtol=1e-10, nans=true)
  y1 = acoustic_field(pm1, tx, rxs; mode=:incoherent)
  y4 = acoustic_field(pm4, tx, rxs; mode=:incoherent)
  @test isapprox(y1, y4; rtol=1e-10, nans=true)
  rx = AcousticReceiver(1000.0, -60.0)
  a1 = arrivals(pm1, tx, rx)
  a4 = arrivals(pm4, tx, rx)
  @test length(a1) == length(a4)
  @test all(isapprox(a1[j].t, a4[j].t; atol=1e-12) for j in eachindex(a1))
  @test all(isapprox(a1[j].ϕ, a4[j].ϕ; rtol=1e-8) for j in eachindex(a1))
end
