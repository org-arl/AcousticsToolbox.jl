using TestItems

@testitem "bellhop-pekeris-rays+ir" begin
  using UnderwaterAcoustics
  env = UnderwaterEnvironment(bathymetry = 20.0u"m", temperature = 27.0u"°C", salinity = 35.0u"ppt", seabed = SandySilt)
  pm = Bellhop(env)
  tx = AcousticSource((x=0.0, z=-5.0), 1000.0)
  rx = AcousticReceiver((x=100.0, z=-10.0))
  arr = arrivals(pm, tx, rx)
  @test arr isa AbstractArray{<:UnderwaterAcoustics.RayArrival}
  @test length(arr) ≥ 6
  @test arr[1].t ≈ 0.0650 atol=0.0001
  @test arr[2].t ≈ 0.0657 atol=0.0001
  @test arr[3].t ≈ 0.0670 atol=0.0001
  @test all([arr[j].t > arr[j-1].t for j ∈ 2:6])
  @test abs(arr[1].ϕ) ≈ 0.01 atol=0.001
  @test real(arr[2].ϕ) < 0.0
  @test imag(arr[2].ϕ) ≈ 0.0 atol=0.001
  @test all([abs(arr[j].ϕ) < abs(arr[j-1].ϕ) for j ∈ 2:6])
  @test [(arr[j].ns, arr[j].nb) for j ∈ 1:6] == [(0,0), (1,0), (0,1), (1,1), (1,1), (2,1)]
  @test all([abs(arr[j].θᵣ) == abs(arr[j].θₛ) for j ∈ 1:6])
  @test all([arr[j].path[1] == (x=0.0, y=0.0, z=-5.0) for j ∈ 1:6])
  @test all([abs(arr[j].path[end].x - 100) < 2 for j ∈ 1:6])
  @test all([arr[j].path[end].y == 0 for j ∈ 1:6])
  @test all([abs(arr[j].path[end].z + 10) < 1 for j ∈ 1:6])
  ir1 = impulse_response(pm, tx, rx, 10000.0; abstime=false)
  ir2 = impulse_response(pm, tx, rx, 10000.0; abstime=true)
  @test length(ir2) ≈ length(ir1) + round(Int, 10000.0 * arr[1].t) atol=1
  @test length(ir2) ≈ round(Int, 10000.0 * arr[end].t) + 1 atol=1
  @test all(abs.(ir1[[round(Int, 10000.0 * (arr[j].t - arr[1].t)) + 1 for j ∈ 1:6]]) .> 0)
  @test all(abs.(ir2[[round(Int, 10000.0 * arr[j].t) + 1 for j ∈ 1:7]]) .> 0)
  @test length(impulse_response(pm, tx, rx, 10000.0; ntaps=256, abstime=true)) == 256
  @test length(impulse_response(pm, tx, rx, 10000.0; ntaps=64, abstime=true)) == 64
  @test length(impulse_response(pm, tx, rx, 10000.0; ntaps=1024, abstime=false)) == 1024
  @test length(impulse_response(pm, tx, rx, 10000.0; ntaps=700, abstime=false)) == 700
end

@testitem "bellhop-pekeris-field+tl" begin
  using UnderwaterAcoustics
  env = UnderwaterEnvironment(bathymetry = 20.0u"m", soundspeed = 1500.0u"m/s",
    seabed = FluidBoundary(water_density(), 1500.0u"m/s", 0.0))  # non-reflecting
  pm = Bellhop(env)
  d = (√1209.0)/4.0
  x = acoustic_field(pm, AcousticSource((x=0.0, z=-d), 1000.0), AcousticReceiver((x=100.0, z=-d)))
  @test x isa Complex
  @test abs(x) ≈ 0.0 atol=0.0002
  x′ = acoustic_field(pm, AcousticSource((x=0.0, z=-d), 1000.0), AcousticReceiver((x=100.0, z=-d)); mode=:incoherent)
  @test x′ isa Complex
  @test imag(x′) == 0.0
  @test abs(x′) > 1/100.0
  d = (√2409.0)/8.0
  x = acoustic_field(pm, AcousticSource((x=0.0, z=-d), 1000.0), AcousticReceiver((x=100.0, z=-d)))
  @test abs(x) > abs(x′)
  y = transmission_loss(pm, AcousticSource((x=0.0, z=-d), 1000.0), AcousticReceiver((x=100.0, z=-d)))
  @test -10 * log10(abs2(x)) ≈ y atol=0.1
  x = acoustic_field(pm, AcousticSource((x=0.0, z=-d), 1000.0), AcousticReceiver((x=100.0, z=-d)); mode=:incoherent)
  @test abs(x) ≈ abs(x′) atol=0.0001
  y = transmission_loss(pm, AcousticSource((x=0.0, z=-d), 1000.0), AcousticReceiver((x=100.0, z=-d)); mode=:incoherent)
  @test -10 * log10(abs2(x)) ≈ y atol=0.1
  tx = AcousticSource((x=0.0, z=-5.0), 1000.0)
  x1 = acoustic_field(pm, tx, AcousticReceiver((x=100.0, z=-5.0)))
  x2 = acoustic_field(pm, tx, AcousticReceiver((x=100.0, z=-10.0)))
  x3 = acoustic_field(pm, tx, AcousticReceiver((x=100.0, z=-15.0)))
  x = acoustic_field(pm, tx, [AcousticReceiver((x=100.0, z=-d)) for d ∈ 5.0:5.0:15.0])
  @test x isa AbstractVector
  @test [x1, x2, x3] == x
  x = acoustic_field(pm, tx, AcousticReceiverGrid2D(100.0:100.0, -5.0:-5.0:-15.0))
  @test x isa AbstractMatrix
  @test size(x) == (1, 3)
  @test [x1 x2 x3] == x
  x = acoustic_field(pm, tx, AcousticReceiverGrid2D(100.0:10.0:120.0, -5.0:-5.0:-15.0))
  @test x isa AbstractMatrix
  @test size(x) == (3, 3)
  @test [x1, x2, x3] == x[1,:]
  x1 = transmission_loss(pm, tx, AcousticReceiver(100.0, -5.0))
  x2 = transmission_loss(pm, tx, AcousticReceiver(100.0, -10.0))
  x3 = transmission_loss(pm, tx, AcousticReceiver(100.0, -15.0))
  x = transmission_loss(pm, tx, [AcousticReceiver(100.0, -d) for d ∈ 5.0:5.0:15.0])
  @test x isa AbstractVector
  @test [x1, x2, x3] == x
  x = transmission_loss(pm, tx, AcousticReceiverGrid2D(100.0:100.0, -5.0:-5.0:-15.0))
  @test x isa AbstractMatrix
  @test size(x) == (1, 3)
  @test [x1 x2 x3] == x
  x = transmission_loss(pm, tx, AcousticReceiverGrid2D(100.0:10.0:120.0, -5.0:-5.0:-15.0))
  @test x isa AbstractMatrix
  @test size(x) == (3, 3)
  @test [x1, x2, x3] == x[1,:]
end

@testitem "bellhop-munk-rays+tl" begin
  using UnderwaterAcoustics
  using Statistics
  env = UnderwaterEnvironment(
    bathymetry = 5.0u"km",
    soundspeed = SampledField(
      [1548.52, 1530.29, 1526.69, 1517.78, 1509.49, 1504.30, 1501.38, 1500.14, 1500.12, 1501.02,
       1502.57, 1504.62, 1507.02, 1509.69, 1512.55, 1515.56, 1518.67, 1521.85, 1525.10, 1528.38,
       1531.70, 1535.04, 1538.39, 1541.76, 1545.14, 1548.52, 1551.91]; z=0:-200:-5000),
    seabed = FluidBoundary(water_density(), 1600.0u"m/s", 0.0))
  pm = Bellhop(env; min_angle=-25°, max_angle=25°, nbeams=5001)
  tx = AcousticSource((x=0, z=-1000u"m"), 50.0)
  rx = AcousticReceiver((x=100u"km", z=-800u"m"))
  arr = arrivals(pm, tx, rx)
  @test length(arr) ≥ 10
  @test arr[1].θₛ ≈ deg2rad(-7.6) atol=0.1
  @test arr[1].ns == 0
  @test arr[1].nb == 0
  @test arr[1].θᵣ ≈ deg2rad(-5.9) atol=0.1
  @test arr[1].t ≈ 66641.04 / 1000 atol = 0.0001
  @test arr[1].ϕ ≈ 10^(-92.8/20) * cispi(180.0/180) atol=0.001
  @test arr[2].θₛ ≈ deg2rad(-5.5) atol=0.1
  @test arr[2].ns == 0
  @test arr[2].nb == 0
  @test arr[2].θᵣ ≈ deg2rad(3.0) atol=0.1
  @test arr[2].t ≈ 66648.70 / 1000 atol = 0.0001
  @test arr[2].ϕ ≈ 10^(-89.7/20) * cispi(-90.0/180) atol=0.001
  @test arr[3].θₛ ≈ deg2rad(-14.2) atol=0.1
  @test arr[3].ns == 1
  @test arr[3].nb == 2
  @test arr[3].θᵣ ≈ deg2rad(-13.5) atol=0.1
  @test arr[3].t ≈ 66825.35 / 1000 atol = 0.0001
  @test arr[3].ϕ ≈ 10^(-105.2/20) * cispi(120.3/180) atol=0.001
  @test arr[4].θₛ ≈ deg2rad(-14.7) atol=0.1
  @test arr[4].ns == 2
  @test arr[4].nb == 2
  @test arr[4].θᵣ ≈ deg2rad(14.0) atol=0.1
  @test arr[4].t ≈ 67017.36 / 1000 atol = 0.0001
  @test arr[4].ϕ ≈ 10^(-104.0/20) * cispi(-85.9/180) atol=0.001
  @test arr[5].θₛ ≈ deg2rad(14.8) atol=0.1
  @test arr[5].ns == 2
  @test arr[5].nb == 2
  @test arr[5].θᵣ ≈ deg2rad(-13.9) atol=0.1
  @test arr[5].t ≈ 67083.31 / 1000 atol = 0.0001
  @test arr[5].ϕ ≈ 10^(-103.7/20) * cispi(-91.1/180) atol=0.001
  @test arr[6].θₛ ≈ deg2rad(15.4) atol=0.1
  @test arr[6].ns == 3
  @test arr[6].nb == 2
  @test arr[6].θᵣ ≈ deg2rad(14.8) atol=0.1
  @test arr[6].t ≈ 67291.82 / 1000 atol = 0.0001
  @test arr[6].ϕ ≈ 10^(-103.1/20) * cispi(63.6/180) atol=0.001
  @test arr[7].θₛ ≈ deg2rad(-17.9) atol=0.1
  @test arr[7].ns == 2
  @test arr[7].nb == 3
  @test arr[7].θᵣ ≈ deg2rad(-17.0) atol=0.1
  @test arr[7].t ≈ 68371.07 / 1000 atol = 0.0001
  @test arr[7].ϕ ≈ 10^(-102.0/20) * cispi(-140.2/180) atol=0.001
  @test arr[8].θₛ ≈ deg2rad(-18.7) atol=0.1
  @test arr[8].ns == 3
  @test arr[8].nb == 3
  @test arr[8].θᵣ ≈ deg2rad(18.2) atol=0.1
  @test arr[8].t ≈ 68652.33 / 1000 atol = 0.0001
  @test arr[8].ϕ ≈ 10^(-101.9/20) * cispi(-9.2/180) atol=0.001
  @test arr[9].θₛ ≈ deg2rad(18.8) atol=0.1
  @test arr[9].ns == 3
  @test arr[9].nb == 3
  @test arr[9].θᵣ ≈ deg2rad(-18.0) atol=0.1
  @test arr[9].t ≈ 68736.42 / 1000 atol = 0.0001
  @test arr[9].ϕ ≈ 10^(-101.8/20) * cispi(-21.9/180) atol=0.001
  @test arr[10].θₛ ≈ deg2rad(19.6) atol=0.1
  @test arr[10].ns == 4
  @test arr[10].nb == 3
  @test arr[10].θᵣ ≈ deg2rad(19.1) atol=0.1
  @test arr[10].t ≈ 69036.50 / 1000 atol = 0.0001
  @test arr[10].ϕ ≈ 10^(-101.8/20) * cispi(85.6/180) atol=0.001
  pm = Bellhop(env; min_angle=-20.3°, max_angle=20.3°, beam_type=:gaussian)
  rxs = AcousticReceiverGrid2D(0:200:100000, -5000:25:0)
  xloss = transmission_loss(pm, tx, rxs; mode=:coherent)
  @test size(xloss) == (501, 201)
  @test minimum(xloss[201:end,:]) ≈ 73.8 atol=0.1
  @test -10 .* log10.(mean(10 .^ (-xloss[201:end,:] ./ 20))) ≈ 44.5 atol=0.1
  xloss = transmission_loss(pm, tx, rxs; mode=:semicoherent)
  @test size(xloss) == (501, 201)
  @test minimum(xloss[201:end,:]) ≈ 77.9 atol=0.1
  @test -10 .* log10.(mean(10 .^ (-xloss[201:end,:] ./ 20))) ≈ 43.8 atol=0.1
  xloss = transmission_loss(pm, tx, rxs; mode=:incoherent)
  @test size(xloss) == (501, 201)
  @test minimum(xloss[201:end,:]) ≈ 77.9 atol=0.1
  @test -10 .* log10.(mean(10 .^ (-xloss[201:end,:] ./ 20))) ≈ 43.8 atol=0.1
end

@testitem "bellhop-seamount-tl" begin
  using UnderwaterAcoustics
  bathy = SampledField([3000, 3000, 500, 3000, 3000]; x=[0, 10e3, 20e3, 30e3, 100e3])
  ssp = SampledField([
    1476.7, 1476.7, 1476.7, 1476.7, 1476.7, 1476.7, 1476.7, 1476.7, 1476.7, 1472.6, 1468.8,
    1467.2, 1471.6, 1473.6, 1473.6, 1472.7, 1472.2, 1471.6, 1471.6, 1472.0, 1472.7, 1473.1,
    1474.9, 1477.0, 1478.1, 1480.7, 1483.8, 1490.5, 1498.3, 1506.5]; z=[0, -5, -10, -15, -20,
    -25, -30, -35, -38, -50, -70, -100, -140, -160, -170, -200, -215, -250, -300, -370, -450,
    -500, -700, -900, -1000, -1250, -1500, -2000, -2500, -3000])
  env = UnderwaterEnvironment(
    bathymetry = bathy,
    soundspeed = ssp,
    seabed = FluidBoundary(1.5*water_density(), 1550.0, dBperλ(0.5)))
  tx = AcousticSource((x=0, z=-18), 230.0)
  rxs = AcousticReceiverGrid2D(0:100:100000, -3000:15:0)
  pm = Bellhop(env; min_angle=-89°, max_angle=89°, beam_type=:gaussian)
  xloss = transmission_loss(pm, tx, rxs; mode=:coherent)
  @test size(xloss) == (1001, 201)
  @test xloss[200,50] > 150
  @test xloss[50,50] ≈ 67.0 atol=0.1
  @test xloss[100,100] ≈ 76.4 atol=0.1
  @test xloss[500,50] ≈ 103.2 atol=0.1
end
