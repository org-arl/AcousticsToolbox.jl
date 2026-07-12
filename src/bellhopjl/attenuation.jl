# SPDX-License-Identifier: GPL-3.0-or-later

# Attenuation conversions. Ports misc/AttenMod.f90 (CRCI, Franc_Garr).

# dB → Nepers conversion constant used throughout Bellhop (20 log10 e)
const DB_PER_NEPER = 8.6858896

"""
    crci(c, α_dbλ, freq; volume=nothing)

Convert real sound speed + attenuation in dB/wavelength ('W' unit — the option
the AcousticsToolbox wrapper writes) into a complex sound speed with positive
imaginary part. Ports `CRCI` (misc/AttenMod.f90):

    αT [Np/m] = α[dB/λ] · f / (8.6858896 · c)  (+ volume attenuation)
    cimag     = αT · c² / ω
    c̃         = c + i·cimag

`volume` may be a `FrancoisGarrison` (wrapper option 'F') to add volume
attenuation, or `:thorp` for Thorp.
"""
function crci(c, α_dbλ, freq; volume=nothing)
    ω = 2π * freq
    αT = iszero(c) ? zero(α_dbλ * freq / c) : α_dbλ * freq / (DB_PER_NEPER * c)
    if volume === :thorp
        αT += thorp(freq) / (DB_PER_NEPER * 1000)              # dB/km → Np/m
    elseif volume isa FrancoisGarrison
        αT += franc_garr(volume, freq / 1000) / (DB_PER_NEPER * 1000)
    end
    cimag = αT * c^2 / ω
    complex(c, cimag)
end

"""
    thorp(freq) -> α in dB/km

Thorp volume attenuation, updated formula from JKPS Eq. 1.34 as in
AttenMod.f90 ('T' option). `freq` in Hz.
"""
function thorp(freq)
    f2 = (freq / 1000)^2
    3.3e-3 + 0.11 * f2 / (1 + f2) + 44 * f2 / (4100 + f2) + 3.0e-4 * f2
end

"Francois-Garrison volume-attenuation parameters (misc/AttenMod.f90 module vars)."
struct FrancoisGarrison{T<:Real}
    temperature::T   # deg C
    salinity::T      # psu
    pH::T
    z_bar::T         # depth [m]
end

FrancoisGarrison(t, s, pH, z̄) = FrancoisGarrison(promote(t, s, pH, z̄)...)

"""
    franc_garr(p, f) -> α in dB/km

Francois-Garrison attenuation at frequency `f` [kHz]. Ports `Franc_Garr`
(misc/AttenMod.f90).
"""
function franc_garr(p::FrancoisGarrison, f)
    T, S, pH, z_bar = p.temperature, p.salinity, p.pH, p.z_bar
    c = 1412 + 3.21 * T + 1.19 * S + 0.0167 * z_bar
    # boric acid contribution
    A1 = 8.86 / c * 10^(0.78 * pH - 5)
    P1 = 1
    f1 = 2.8 * sqrt(S / 35) * 10^(4 - 1245 / (T + 273))
    # magnesium sulfate contribution
    A2 = 21.44 * S / c * (1 + 0.025 * T)
    P2 = 1 - 1.37e-4 * z_bar + 6.2e-9 * z_bar^2
    f2 = 8.17 * 10^(8 - 1990 / (T + 273)) / (1 + 0.0018 * (S - 35))
    # viscosity
    P3 = 1 - 3.83e-5 * z_bar + 4.9e-10 * z_bar^2
    A3 = T < 20 ?
        4.937e-4 - 2.59e-5 * T + 9.11e-7 * T^2 - 1.5e-8 * T^3 :
        3.964e-4 - 1.146e-5 * T + 1.45e-7 * T^2 - 6.5e-10 * T^3
    A1 * P1 * (f1 * f^2) / (f1^2 + f^2) + A2 * P2 * (f2 * f^2) / (f2^2 + f^2) +
        A3 * P3 * f^2
end
