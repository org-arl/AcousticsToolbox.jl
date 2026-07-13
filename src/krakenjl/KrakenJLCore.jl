# SPDX-License-Identifier: GPL-3.0-or-later
#
# KrakenJLCore — native Julia port of the KRAKEN / KRAKENC normal-mode models
# (derivative work of KRAKEN by Michael B. Porter, from the OALIB Acoustics
# Toolbox release 2024_12_25). Unlike the rest of AcousticsToolbox.jl (MIT),
# everything under src/krakenjl/ is licensed GPL-3.0-or-later — see
# LICENSE-GPL-3.0 at the repository root.
#
# This internal submodule keeps the port's helpers out of the parent
# namespace; the public model type `KrakenJL` is re-exported by
# AcousticsToolbox. Comments cite the originating Fortran file/subroutine
# throughout so the code remains diffable against upstream. See
# PORTING_NOTES.md in this directory for lineage and deviations.

module KrakenJLCore

using LinearAlgebra

include("attenuation.jl")
include("types.jl")
include("dual.jl")
include("mesh.jl")
include("bcimp.jl")
include("solve.jl")
include("modes.jl")
include("field.jl")
include("api.jl")

export KrakenJL

end # module
