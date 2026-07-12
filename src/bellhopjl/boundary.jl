# SPDX-License-Identifier: GPL-3.0-or-later

# Piecewise-linear altimetry/bathymetry. Ports bdryMod.f90
# (ComputeBdryTangentNormal, GetTopSeg/GetBotSeg). Normals are OUTWARD
# pointing (out of the water column), as in the Fortran:
#   bottom: n = (-t₂, +t₁)  (flat bottom ⇒ (0, +1), pointing down/out)
#   top:    n = (+t₂, -t₁)  (flat top    ⇒ (0, -1), pointing up/out)
# Distances2D: Dist = -dot(n, x - bdry_x) is positive inside the water.

struct Boundary2D{T<:Real}
    rnode::Vector{T}          # node ranges (increasing, extended to ±BDRY_BIG)
    x::Vector{SVector{2,T}}   # nodes (r, z)
    t::Vector{SVector{2,T}}   # unit tangent per segment (length n-1)
    n::Vector{SVector{2,T}}   # unit outward normal per segment
    len::Vector{T}            # segment lengths
    κ::Vector{T}              # curvature per segment (0 for piecewise linear)
end

"""
    Boundary2D(r, z; top)

Piecewise-linear boundary through nodes `(r[i], z[i])`, extended to ±∞ in a
piecewise-constant fashion as in `ComputeBdryTangentNormal` (bdryMod.f90).
"""
function Boundary2D(r::AbstractVector, z::AbstractVector; top::Bool)
    T = float(promote_type(eltype(r), eltype(z)))
    m = length(r)
    @assert length(z) == m ≥ 1
    big = T(BDRY_BIG)
    rv = vcat(-big, collect(T, r), big)
    zv = vcat(T(z[1]), collect(T, z), T(z[end]))
    n = length(rv)
    pts = [SVector{2,T}(rv[i], zv[i]) for i in 1:n]
    tv = Vector{SVector{2,T}}(undef, n - 1)
    nv = Vector{SVector{2,T}}(undef, n - 1)
    lv = Vector{T}(undef, n - 1)
    for i in 1:n-1
        d = pts[i+1] - pts[i]
        l = norm(d)
        d = d / l
        tv[i] = d
        lv[i] = l
        nv[i] = top ? SVector(d[2], -d[1]) : SVector(-d[2], d[1])
    end
    Boundary2D{T}(rv, pts, tv, nv, lv, zeros(T, n - 1))
end

"Flat boundary at depth z0."
flat_boundary(z0; top::Bool) = Boundary2D([zero(z0)], [z0]; top)

"""
    get_seg(b, r, tr) -> i

Segment index for range `r` with the tangent-direction-aware edge handling of
`GetTopSeg`/`GetBotSeg` (bdryMod.f90): with `tr > 0` the segment satisfies
`rnode[i] <= r`, otherwise `rnode[i] < r`.
"""
@inline function get_seg(b::Boundary2D, r, tr)
    n = length(b.rnode)
    i = tr > 0 ? searchsortedlast(b.rnode, r) :
                 searchsortedfirst(b.rnode, r) - 1   # last node with rnode < r
    clamp(i, 1, n - 1)
end

"Geometry of segment i: anchor point, tangent, outward normal, r-interval, curvature."
@inline function segment_geom(b::Boundary2D, i::Int)
    (x = b.x[i], t = b.t[i], n = b.n[i], len = b.len[i],
     rseg = (b.rnode[i], b.rnode[i+1]), κ = b.κ[i])
end

"""
Signed perpendicular distances of x to top and bottom (positive = inside the
water). Ports `Distances2D` (bellhop.f90) — outward normals, hence the minus.
"""
@inline function distances(x::SVector{2}, topg, botg)
    dtop = -dot(topg.n, x - topg.x)
    dbot = -dot(botg.n, x - botg.x)
    dtop, dbot
end
