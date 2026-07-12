# SPDX-License-Identifier: GPL-3.0-or-later

"""
Internal, self-contained 2D environment in Bellhop-native conventions
(z positive down, source at r = 0). Constructed by the API adapter (api.jl).

The numeric element types of the components may differ (e.g. dual-number SSP
nodes with plain-Float64 boundaries when differentiating w.r.t. soundspeed);
`env_eltype` gives the promoted scalar type.
"""
struct Env2D{T1<:Real, S<:AbstractSSP, B1<:Boundary2D, B2<:Boundary2D,
             TB<:BoundaryCondition, BB<:BoundaryCondition}
    freq::T1
    ssp::S
    top::B1
    bot::B2
    topbc::TB
    botbc::BB
end

Base.eltype(::Boundary2D{T}) where {T} = T
bc_eltype(::BoundaryCondition) = Bool          # neutral in promote_type
bc_eltype(::HalfspaceBC{T}) where {T} = T
bc_eltype(::TabulatedBC{T}) where {T} = T

env_eltype(env::Env2D) =
    promote_type(typeof(float(env.freq)), eltype(env.ssp), eltype(env.top),
                 eltype(env.bot), bc_eltype(env.topbc), bc_eltype(env.botbc),
                 Float64)

"Evaluate the SSP at position x = (r, z) with scaled tangent t. Returns (SSPEval, iseg)."
@inline evalssp(env::Env2D, x::SVector{2}, t::SVector{2}, iseg::Int) =
    evaluate(env.ssp, x[2], t[2], iseg)
