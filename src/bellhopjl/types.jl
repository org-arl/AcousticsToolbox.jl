# SPDX-License-Identifier: GPL-3.0-or-later

"""
    RayPt{T}

State of one point along a ray. Mirrors `ray2DPt` in `bellhopMod.f90`.

- `x` — position (r, z), z positive down [m]
- `t` — *scaled* tangent, `t = unit_tangent / c` (so `|t| = 1/c`) [s/m]
- `p, q` — dynamic ray-tracing quantities; initial conditions p=(1,0), q=(0,1)
  (geometric beams use only the first component and start with q=(0,0))
- `τ` — complex travel time; imaginary part accumulates volume attenuation
- `c` — sound speed at `x` [m/s]
- `amp`, `phase` — ray amplitude and accumulated phase from boundary reflections
- `ntop`, `nbot` — boundary bounce counters
"""
struct RayPt{T<:Real}
    x::SVector{2,T}
    t::SVector{2,T}
    p::SVector{2,T}
    q::SVector{2,T}
    τ::Complex{T}
    c::T
    amp::T
    phase::T
    ntop::Int
    nbot::Int
end

"Beam/run configuration. Mirrors the `Beam` structure in `bellhopMod.f90`."
Base.@kwdef struct BeamParams{T<:Real}
    deltas::T                 # step size [m] (0 not allowed here; resolved by caller)
    box_r::T                  # ray box max range [m]
    box_z::T                  # ray box max depth [m]
    nbeams::Int
    αmin::T                   # take-off angle limits [rad], positive down
    αmax::T
    gaussian::Bool = false    # false ⇒ geometric hat ('G'), true ⇒ Gaussian ('B')
end

@enum RunMode COHERENT INCOHERENT SEMICOHERENT ARRIVALS

# bellhopMod.f90: MaxN = 100000 (max # of steps along a ray);
# TraceRay2D exits at is >= MaxN - 3.
const MAX_N = 100_000

# Step.f90: INFINITESIMAL_STEP_SIZE = 1.0d-6
const INFINITESIMAL_STEP_SIZE = 1.0e-6

# bdryMod.f90 / ComputeBdryTangentNormal: boundaries are extended to
# ±sqrt(huge)/1e5 in range
const BDRY_BIG = sqrt(floatmax(Float64)) / 1.0e5
