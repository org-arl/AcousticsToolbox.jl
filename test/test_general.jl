using TestItems

@testitem "models" begin
  using UnderwaterAcoustics
  m = models()
  @test m isa Vector{Type{<:UnderwaterAcoustics.AbstractPropagationModel}}
  @test Bellhop ∈ m
  @test Kraken ∈ m
  env = UnderwaterEnvironment(bathymetry = 20.0u"m", temperature = 27.0u"°C", salinity = 35.0u"ppt", seabed = SandySilt)
  m = models(env)
  @test m isa Vector{Type{<:UnderwaterAcoustics.AbstractPropagationModel}}
  @test Bellhop ∈ m
  @test Kraken ∈ m
end
