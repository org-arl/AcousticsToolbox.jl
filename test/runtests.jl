using AcousticsToolbox
using AcousticsToolbox.UnderwaterAcoustics
using Test

@testset "pm-bellhop" begin

  env = UnderwaterEnvironment(seasurface=Vacuum)
  pm = Bellhop(env)
  @test pm isa Bellhop

  arr = arrivals(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiver(100.0, -10.0))
  @test arr isa AbstractArray{<:UnderwaterAcoustics.RayArrival}
  r = eigenrays(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiver(100.0, -10.0))
  @test r isa AbstractArray{<:UnderwaterAcoustics.RayArrival}
  r = rays(pm, AcousticSource(0.0, -5.0, 1000.0), -60°:15°:60°, 100.0)
  @test r isa AbstractArray{<:UnderwaterAcoustics.RayArrival}
  x = transfercoef(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiver(100.0, -5.0))
  @test x isa Complex
  y = transmissionloss(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiver(100.0, -5.0))
  @test -10 * log10(abs2(x)) ≈ y atol=0.1
  x = transfercoef(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiver(100.0, -5.0); mode=:incoherent)
  @test x isa Complex
  @test imag(x) == 0.0
  y = transmissionloss(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiver(100.0, -5.0); mode=:incoherent)
  @test -10 * log10(abs2(x)) ≈ y atol=0.1
  x1 = transfercoef(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiver(100.0, -5.0))
  x2 = transfercoef(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiver(100.0, -10.0))
  x3 = transfercoef(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiver(100.0, -15.0))
  x = transfercoef(pm, AcousticSource(0.0, -5.0, 1000.0), [AcousticReceiver(100.0, -d) for d ∈ 5.0:5.0:15.0])
  @test x isa AbstractVector
  @test [x1, x2, x3] == x
  x = transfercoef(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiverGrid2D(100.0, 0.0, 1, -5.0, -5.0, 3))
  @test x isa AbstractMatrix
  @test size(x) == (1, 3)
  @test [x1 x2 x3] == x
  x = transfercoef(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiverGrid2D(100.0, 10.0, 3, -5.0, -5.0, 3))
  @test x isa AbstractMatrix
  @test size(x) == (3, 3)
  @test [x1, x2, x3] == x[1,:]
  x1 = transmissionloss(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiver(100.0, -5.0))
  x2 = transmissionloss(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiver(100.0, -10.0))
  x3 = transmissionloss(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiver(100.0, -15.0))
  x = transmissionloss(pm, AcousticSource(0.0, -5.0, 1000.0), [AcousticReceiver(100.0, -d) for d ∈ 5.0:5.0:15.0])
  @test x isa AbstractVector
  @test [x1, x2, x3] == x
  x = transmissionloss(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiverGrid2D(100.0, 0.0, 1, -5.0, -5.0, 3))
  @test x isa AbstractMatrix
  @test size(x) == (1, 3)
  @test [x1 x2 x3] == x
  x = transmissionloss(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiverGrid2D(100.0, 10.0, 3, -5.0, -5.0, 3))
  @test x isa AbstractMatrix
  @test size(x) == (3, 3)
  @test [x1, x2, x3] == x[1,:]

  tx = AcousticSource(0.0, -5.0, 1000.0)
  rx = AcousticReceiver(100.0, -10.0)
  sig = record(pm, tx, rx, 1.0, 44100.0)
  @test size(sig) == (44100,)
  sig = record(pm, tx, rx, 1.0, 44100.0; start=0.5)
  @test size(sig) == (44100,)
  sig = recorder(pm, tx, rx)(1.0, 44100.0)
  @test size(sig) == (44100,)
  sig = recorder(pm, tx, rx)(1.0, 44100.0; start=0.5)
  @test size(sig) == (44100,)
  tx = [AcousticSource(0.0, -5.0, 1000.0), AcousticSource(0.0, -10.0, 2000.0)]
  rx = AcousticReceiver(100.0, -10.0)
  sig = record(pm, tx, rx, 1.0, 44100.0)
  @test size(sig) == (44100,)
  sig = record(pm, tx, rx, 1.0, 44100.0; start=0.5)
  @test size(sig) == (44100,)
  sig = recorder(pm, tx, rx)(1.0, 44100.0)
  @test size(sig) == (44100,)
  sig = recorder(pm, tx, rx)(1.0, 44100.0; start=0.5)
  @test size(sig) == (44100,)
  tx = [AcousticSource(0.0, -5.0, 1000.0), AcousticSource(0.0, -10.0, 2000.0)]
  rx = [AcousticReceiver(100.0, -10.0), AcousticReceiver(100.0, -15.0)]
  sig = record(pm, tx, rx, 1.0, 44100.0)
  @test size(sig) == (44100,2)
  sig = record(pm, tx, rx, 1.0, 44100.0; start=0.5)
  @test size(sig) == (44100,2)
  sig = recorder(pm, tx, rx)(1.0, 44100.0)
  @test size(sig) == (44100,2)
  sig = recorder(pm, tx, rx)(1.0, 44100.0; start=0.5)
  @test size(sig) == (44100,2)
  tx = AcousticSource(0.0, -5.0, 1000.0)
  rx = [AcousticReceiver(100.0, -10.0), AcousticReceiver(100.0, -15.0)]
  sig = record(pm, tx, rx, 1.0, 44100.0)
  @test size(sig) == (44100,2)
  sig = record(pm, tx, rx, 1.0, 44100.0; start=0.5)
  @test size(sig) == (44100,2)
  sig = recorder(pm, tx, rx)(1.0, 44100.0)
  @test size(sig) == (44100,2)
  sig = recorder(pm, tx, rx)(1.0, 44100.0; start=0.5)
  @test size(sig) == (44100,2)

  env = UnderwaterEnvironment(
    seasurface=Vacuum,
    ssp=SampledSSP(0.0:5.0:20.0, [1500.0, 1490.0, 1500.0, 1505.0, 1507.0]),
    altimetry=SampledAltitude(0.0:25.0:100.0, [0.0, -1.0, 0.0, -1.0, 0.0]),
    bathymetry=SampledDepth(0.0:25.0:100.0, [20.0, 17.0, 17.0, 19.0, 20.0])
  )
  pm = Bellhop(env)
  @test pm isa Bellhop
  r = eigenrays(pm, AcousticSource(0.0, -5.0, 5000.0), AcousticReceiver(100.0, -10.0))
  @test r isa AbstractArray{<:UnderwaterAcoustics.RayArrival}
  x = transmissionloss(pm, AcousticSource(0.0, -5.0, 5000.0), AcousticReceiverGrid2D(1.0, 1.0, 100, 0.0, -1.0, 20))
  @test x isa AbstractMatrix
  @test size(x) == (100, 20)

  struct TestAlt <: Altimetry end
  UnderwaterAcoustics.altitude(::TestAlt, x, y) = -1.0 + sin(2π*x/10.0)

  struct TestBathy <: Bathymetry end
  UnderwaterAcoustics.depth(::TestBathy, x, y) = 18.0 + 2*sin(2π*x/30.0)
  UnderwaterAcoustics.maxdepth(::TestBathy) = 20.0

  env = UnderwaterEnvironment(
    seasurface=Vacuum,
    ssp=MunkSSP(),
    altimetry=TestAlt(),
    bathymetry=TestBathy()
  )
  pm = Bellhop(env)
  @test pm isa Bellhop
  r = eigenrays(pm, AcousticSource(0.0, -5.0, 5000.0), AcousticReceiver(100.0, -10.0))
  @test r isa AbstractArray{<:UnderwaterAcoustics.RayArrival}
  x = transmissionloss(pm, AcousticSource(0.0, -5.0, 5000.0), AcousticReceiverGrid2D(1.0, 1.0, 100, 0.0, -1.0, 20))
  @test x isa AbstractMatrix
  @test size(x) == (100, 20)

end
