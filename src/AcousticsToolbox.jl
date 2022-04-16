module AcousticsToolbox

import AcousticsToolbox_jll

using UnderwaterAcoustics
using UnderwaterAcoustics: RayArrival, SampledSSP1D, check
using Printf

include("Bellhop.jl")
include("Kraken.jl")

function __init__()
  UnderwaterAcoustics.addmodel!(Bellhop)
  UnderwaterAcoustics.addmodel!(Kraken)
end

end # module
