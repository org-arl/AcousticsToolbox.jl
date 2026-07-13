# SPDX-License-Identifier: GPL-3.0-or-later

# Problem description and solver workspace. Mirrors the module variables of
# KrakenMod.f90 / KrakencMod.f90 and the SSP structure of misc/sspMod.f90,
# repackaged as explicit structs (no globals).

# KrakenMod.f90: MaxM = 20000 (max # of modes), NSets = 5 (Richardson meshes)
const MAX_MODES = 20000
const NSETS = 5
const NV = (1, 2, 4, 8, 16)          # mesh multipliers, kraken(c).f90 main

# FUNCT / AcousticLayers / BCImpedance scaling constants (kraken(c).f90)
const IPOW_R = 50
const IPOW_F = -50
const ROOF = 1.0e50
const FLOOR_ = 1.0e-50

"""
Halfspace / boundary condition. Mirrors `HSInfo` (sspMod.f90) restricted to
the codes reachable through the wrapper: 'V' vacuum, 'R' rigid, 'A'
acousto-elastic halfspace (fluid when cs == 0).
"""
struct Halfspace{T<:Real}
    bc::Symbol               # :vacuum, :rigid, :acoustoelastic
    cp::Complex{T}           # complex compressional speed (CRCI-converted)
    cs::Complex{T}           # complex shear speed (0 for fluid)
    rho::T                   # density ratio (g/cm³, water = 1)
end

Halfspace{T}(bc::Symbol) where {T} =
    Halfspace{T}(bc, zero(Complex{T}), zero(Complex{T}), zero(T))

"""
One medium (layer) of the problem. SSP nodes are CRCI-converted complex
speeds; `spline` selects cubic-spline ('S') vs C-linear ('C') subtabulation
(sspMod.f90 cCubic/cLinear). `ng` is the base finite-difference mesh count
(NG in the Fortran; the auto rule of ReadEnvironmentMod is applied by the
caller). `sigma` roughness values live in `Problem` (per interface).
"""
struct Medium{T<:Real}
    z0::T                    # depth of top of medium [m]
    z1::T                    # depth of bottom of medium [m]
    zn::Vector{T}            # SSP node depths (ascending, zn[1]=z0, zn[end]=z1)
    cpn::Vector{Complex{T}}  # complex P-speed at nodes
    csn::Vector{Complex{T}}  # complex S-speed at nodes (nonzero ⇒ elastic)
    rhon::Vector{T}          # density at nodes (g/cm³ ratio)
    spline::Bool
    ng::Int
end

is_elastic(m::Medium) = any(!iszero, m.csn)

"Range-independent normal-mode problem (one profile, one frequency)."
struct Problem{T<:Real}
    freq::T
    media::Vector{Medium{T}}
    sigma::Vector{T}         # interface roughness, length nmedia+1 (top..bottom)
    hstop::Halfspace{T}
    hsbot::Halfspace{T}
    clow::T                  # phase-speed search interval
    chigh::T
    rmax::T                  # max range [m] for the Richardson convergence test
end

"""
Finite-difference workspace for one mesh (one `iSet`). Mirrors the module
arrays of Kraken(c)Mod.f90. `X` is the arithmetic type of the characteristic
function: `Complex{T}` for KRAKENC, `T` for KRAKEN (which keeps the imaginary
part of ω²/c² separately in `b1c`).
"""
mutable struct MeshWorkspace{T<:Real,X<:Number}
    n::Vector{Int}           # mesh intervals per medium
    h::Vector{T}             # mesh spacing per medium
    loc::Vector{Int}         # offset of each medium in the point arrays
    b1::Vector{X}            # FD diagonal (acoustic) / elastic coefficient
    b1c::Vector{T}           # KRAKEN only: imag(ω²/c²) for the perturbation
    b2::Vector{X}            # elastic coefficients (kraken(c).f90 Initialize)
    b3::Vector{X}
    b4::Vector{X}
    rho::Vector{T}
    first_acoustic::Int
    last_acoustic::Int
    cmin::T
    clow::T                  # search interval after the cmin adjustments
    chigh::T
end

function MeshWorkspace{T,X}(nmedia::Int) where {T,X}
    MeshWorkspace{T,X}(zeros(Int, nmedia), zeros(T, nmedia), zeros(Int, nmedia),
                       X[], T[], X[], X[], X[], T[], 0, 0, zero(T), zero(T),
                       zero(T))
end

"Result of the mode solve: everything the API layer needs."
struct ModeResult{T<:Real}
    k::Vector{Complex{T}}    # modal wavenumbers (perturbed, sorted)
    phi::Matrix{Complex{T}}  # mode shapes on the solver mesh (NTotal1 × M)
    z::Vector{T}             # mesh depths for phi (acoustic media only)
    vg::Vector{T}            # group speeds (0 when invalid — KRAKEN real path)
end
