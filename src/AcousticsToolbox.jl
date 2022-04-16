module AcousticsToolbox

import AcousticsToolbox_jll

using UnderwaterAcoustics
using UnderwaterAcoustics: RayArrival, SampledSSP1D, check
using Printf

include("Bellhop.jl")

function __init__()
  UnderwaterAcoustics.addmodel!(Bellhop)
end

end # module
