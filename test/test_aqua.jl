using TestItems

@testitem "aqua" begin
  using AcousticsToolbox
  using Aqua
  Aqua.test_all(AcousticsToolbox; persistent_tasks=false)
end
