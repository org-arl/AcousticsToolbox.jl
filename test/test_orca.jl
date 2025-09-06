using TestItems

@testitem "orca-pekeris-modes" begin
  using UnderwaterAcoustics
  env = @inferred UnderwaterEnvironment(
    bathymetry = 100,
    soundspeed = 1500,
    density = 1000,
    seabed = FluidBoundary(2000, 2000)
  )
  tx = @inferred AcousticSource(0, -50, 300)
  rx = @inferred AcousticReceiver(20000, -25)
  for complex_solver ∈ [false, true]
    pm = @inferred Orca(env; complex_solver)
    m = @inferred arrivals(pm, tx, rx)
    @test m isa Vector{<:UnderwaterAcoustics.ModeArrival}
    pm1 = @inferred PekerisModeSolver(env)
    m1 = @inferred arrivals(pm1, tx, rx)
    @test [m[i].kᵣ for i ∈ 1:25] ≈ [m1[i].kᵣ for i ∈ 1:25] atol=0.001
    @test [m[i].v for i ∈ 1:25] ≈ [m1[i].v for i ∈ 1:25] atol=0.1
    rxs = @inferred AcousticReceiverGrid2D(20000:10:22000, -25)
    x = @inferred transmission_loss(pm, tx, rxs)
    x1 = @inferred transmission_loss(pm1, tx, rxs)
    @test sum(abs, x1 - x) / length(x) < 1
  end
end
