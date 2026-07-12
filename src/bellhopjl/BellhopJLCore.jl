# SPDX-License-Identifier: GPL-3.0-or-later
#
# BellhopJLCore — native Julia port of the 2D Bellhop Gaussian-beam/ray tracer
# (derivative work of Bellhop by Michael B. Porter, via the A-New-BellHope
# bug-fixed mirror of Bellhop 2022_4). Unlike the rest of AcousticsToolbox.jl
# (MIT), everything under src/bellhopjl/ is licensed GPL-3.0-or-later — see
# LICENSE-GPL-3.0 at the repository root.
#
# This internal submodule keeps the port's helpers out of the parent
# namespace; the public model type `BellhopJL` is re-exported by
# AcousticsToolbox.

module BellhopJLCore

using LinearAlgebra
using StaticArrays

include("types.jl")
include("ssp.jl")
include("boundary.jl")
include("attenuation.jl")
include("reflection.jl")
include("env.jl")
include("step.jl")
include("trace.jl")
include("influence.jl")
include("field.jl")
include("api.jl")

export BellhopJL

end # module
