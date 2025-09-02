using TestItems

@testitem "models" begin
  using UnderwaterAcoustics
  m = models()
  @test m isa Vector{Type{<:UnderwaterAcoustics.AbstractPropagationModel}}
  @test Bellhop ∈ m
  @test Kraken ∈ m
end
