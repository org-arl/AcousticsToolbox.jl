using TestItems

@testitem "kraken-pekeris-modes" begin
  using UnderwaterAcoustics
  env = UnderwaterEnvironment(
    bathymetry = 5000,
    soundspeed = 1500,
    density = 1000,
    seabed = FluidBoundary(2000, 2000)
  )
  for leaky ∈ [false, true]
    pm = @inferred Kraken(env; chigh=2000, leaky)
    tx = AcousticSource(0, -500, 10)
    rx = AcousticReceiver(200000, -2500)
    m = @inferred arrivals(pm, tx, rx)
    @test m isa Vector{<:UnderwaterAcoustics.ModeArrival}
    pm1 = PekerisModeSolver(env)
    m1 = arrivals(pm1, tx, rx)
    @test length(m) == length(m1)
    @test all(m[i].kᵣ ≈ m1[i].kᵣ for i ∈ eachindex(m))
    # we don't test group velocity because kraken seems to give all 0
    rxs = AcousticReceiverGrid2D(200000:10:220000, -2500)
    x = @inferred transmission_loss(pm, tx, rxs)
    x1 = transmission_loss(pm1, tx, rxs)
    @test sum(abs, x1 - x) / length(x) < 1
  end
end
