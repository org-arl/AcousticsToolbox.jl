module AcousticsToolbox

import AcousticsToolbox_jll

using UnderwaterAcoustics
using UnderwaterAcoustics: AbstractRayPropagationModel, AbstractModePropagationModel
using UnderwaterAcoustics: AbstractAcousticSource, AbstractAcousticReceiver, AcousticReceiverGrid2D
using UnderwaterAcoustics: RayArrival, ModeArrival, SampledFieldX, SampledFieldXZ, SampledFieldZ
using UnderwaterAcoustics: is_range_dependent, is_constant, value, in_units, db2amp
using Printf

include("bellhop.jl")
include("kraken.jl")
include("orca.jl")
include("common.jl")

# BellhopJL (native Julia port of Bellhop) lives in an internal submodule;
# note that unlike the rest of this package (MIT), src/bellhopjl/** is
# GPL-3.0-or-later (see LICENSE-GPL-3.0)
include("bellhopjl/BellhopJLCore.jl")
using .BellhopJLCore: BellhopJL
export BellhopJL

end # module
