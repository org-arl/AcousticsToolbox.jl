# SPDX-License-Identifier: GPL-3.0-or-later

# Boundary conditions and reflection coefficients. Ports the boundary-condition
# branch of Reflect2D (bellhop.f90) and InterpolateReflectionCoefficient
# (misc/RefCoef.f90).

abstract type BoundaryCondition end

struct VacuumBC <: BoundaryCondition end                            # 'V'
struct RigidBC  <: BoundaryCondition end                            # 'R'

"Fluid halfspace ('A'). `ρ` is the density *ratio* (g/cm³ with water = 1)."
struct HalfspaceBC{T<:Real} <: BoundaryCondition
    ρ::T           # density ratio (halfspace ρ / water ρ)
    cp::Complex{T} # complex compressional speed (CRCI of c, α[dB/λ])
end

function HalfspaceBC(ρ, c, α_dbλ, freq)
    cp = crci(c, α_dbλ, freq)
    T = promote_type(typeof(float(ρ)), typeof(real(cp)))
    HalfspaceBC{T}(T(ρ), Complex{T}(cp))
end

"Tabulated reflection coefficient ('F'): grazing angle [deg], |R|, phase [rad]."
struct TabulatedBC{T<:Real} <: BoundaryCondition
    θ::Vector{T}
    R::Vector{T}
    φ::Vector{T}
end

"""
    reflection_coefficient(bc, Tg, Th, ω) -> (amp_scale, phase_add)

Reflection amplitude/phase updates for a boundary, given the tangential (`Tg`)
and normal (`Th`) components of the *scaled* incident tangent along the
boundary tangent/outward normal. Ports the `SELECT CASE ( HS%BC )` block of
`Reflect2D` (bellhop.f90). Returns the multiplicative amplitude factor and
additive phase.
"""
reflection_coefficient(::RigidBC, Tg, Th, ω) = (one(Tg), zero(Tg))
reflection_coefficient(::VacuumBC, Tg, Th, ω) = (one(Tg), oftype(Tg, π))

function reflection_coefficient(bc::HalfspaceBC, Tg, Th, ω)
    kx = ω * Tg           # wavenumber along boundary
    kz = ω * Th           # wavenumber perpendicular to boundary (in ocean)
    # fluid bottom branch (HS%cS = 0): kzP = sqrt(kx² − (ω/cP)²)
    kzP = sqrt(complex(kx^2) - (ω / bc.cp)^2)
    # Intel/GFortran differ on the branch of sqrt for negative reals:
    if real(kzP) == 0 && imag(kzP) < 0
        kzP = -kzP
    end
    f = kzP
    g = bc.ρ
    ρw = one(real(f))     # water density is 1 in Bellhop's units
    R = -(ρw * f - im * kz * g) / (ρw * f + im * kz * g)
    aR = abs(R)
    if aR < 1.0e-5        # kill a ray that has lost its energy in reflection
        (zero(aR), zero(aR))
    else
        (aR, atan(imag(R), real(R)))
    end
end

function reflection_coefficient(bc::TabulatedBC, Tg, Th, ω)
    # angle of incidence relative to normal; symmetric about 90°
    θ = abs(atand(Th, Tg))
    θ > 90 && (θ = 180 - θ)
    # InterpolateReflectionCoefficient (misc/RefCoef.f90): linear interpolation,
    # R = 0 outside the tabulated domain
    n = length(bc.θ)
    i = clamp(searchsortedlast(bc.θ, θ), 1, n - 1)
    if θ < bc.θ[1] || θ > bc.θ[n]
        return (zero(θ), zero(θ))
    end
    w = (θ - bc.θ[i]) / (bc.θ[i+1] - bc.θ[i])
    Rm = bc.R[i] + w * (bc.R[i+1] - bc.R[i])
    φ = bc.φ[i] + w * (bc.φ[i+1] - bc.φ[i])
    (Rm, φ)
end
