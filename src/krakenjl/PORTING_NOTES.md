# KrakenJL porting notes — provenance and deviations

This directory contains a native Julia port of the KRAKEN and KRAKENC normal
mode models. This file documents the port's lineage, every intentional
deviation from the Fortran sources, and the provenance of the validation
data. (Everything in `src/krakenjl/` is GPL-3.0-or-later — see
`LICENSE-GPL-3.0` at the repository root; the rest of AcousticsToolbox.jl is
MIT.)

## 1. Lineage

Ported from the **OALIB Acoustics Toolbox release 2024_12_25**
(`http://oalib.hlsresearch.com/AcousticsToolbox/at_2024_12_25.zip`, sha256
`7b57e80bded7f71ea9536e541029615f3f430e390651d697a2212569cbafd85c`) — the
exact source tree the `AcousticsToolbox_jll` 2025.9.6 binaries are compiled
from (per the Yggdrasil build recipe), so the port and the wrapped Fortran
`Kraken` model share their algorithm version. KRAKEN © Michael B. Porter,
GPL-3.0.

Files ported (docstrings cite the originating subroutine throughout):

| Julia file | Fortran source |
|---|---|
| `types.jl` | `KrakenMod.f90`, `KrakencMod.f90` (module variables → structs) |
| `attenuation.jl` | `misc/AttenMod.f90` (CRCI, Franc_Garr) |
| `mesh.jl` | `Initialize` (kraken(c).f90), `cLinear`/`cCubic` (misc/sspMod.f90), `CSPLINE` (misc/splinec.f90), auto-mesh rule of `ReadEnvironmentMod.f90` |
| `bcimp.jl` | `BCImpedanceMod.f90`, `BCImpedancecMod.f90`, `misc/PekRoot.f90` |
| `solve.jl` | `FUNCT`/`AcousticLayers`/`Solve`/`Solve1`/`Solve2`/`Bisection` (kraken(c).f90), `RootFinderBrent.f90`, `misc/RootFinderSecantMod.f90` |
| `modes.jl` | `Vector`/`Normalize`/`ScatterLoss` (kraken(c).f90), `InverseIterationMod.f90`, `KupIng` (Kraken/Scattering.f90) |
| `field.jl` | `Evaluate` (KrakenField/EvaluateMod.f90), `Weight` (misc/calculateweights.f90) |
| `api.jl` | behavioral port of the wrapper `src/kraken.jl` + `src/common.jl _write_env` |

The ORCA sources were **not** consulted (restrictive ARL:UT license).

## 2. Scope

Only the features reachable through the UnderwaterAcoustics.jl API of the
Fortran wrapper are implemented: one profile (range-independent), one
frequency, point source, cylindrical coordinates, coherent/incoherent mode
summation ('R… C/I' field options), fluid/elastic/multilayer-elastic
seabeds, vacuum/rigid/acousto-elastic surface, C-linear and cubic-spline
SSPs, 'WF' attenuation (dB/λ + Francois-Garrison volume attenuation). Not
ported: tabulated/precalculated reflection coefficients ('F'/'P' BCs),
broadband mode files, profile continuation (`Solve3`), BOUNCE, field3d,
coupled modes, twersky scatter, biological attenuation.

## 3. Intentional deviations from the Fortran

- **No env/mod/shd/flp file I/O** — environments come directly from
  UnderwaterAcoustics.jl; the adapter (`api.jl`) mirrors the unit conversions
  of this package's Fortran env-file writer verbatim (density ratio in g/cm³,
  attenuation in dB/λ, Francois-Garrison volume attenuation with
  z̄ = waterdepth/2 applied through CRCI to every medium and halfspace, mesh
  count `round(waterdepth/λ·mesh_density)` with the ReadEnvironmentMod
  auto rule when 0, `rmax = Inf → 1.01·max receiver range`, converted to km
  for the Richardson convergence test).
- **Precision** — double precision throughout. The Fortran writes mode
  shapes, eigenvalues and the pressure field in single precision (`.mod` /
  `.shd` records are `COMPLEX*4`) and the field program sums modes in single
  precision; this port keeps `Float64`, so results differ from the wrapped
  binaries at the ~1e-7 relative level.
- **Mode shapes on the solver mesh** — `arrivals` returns ψ sampled on the
  solver's own finite-difference mesh instead of the wrapper's λ/10
  resampling (finer; no information loss). Field evaluation tabulates ψ at
  source/receiver depths by the same linear interpolation as `Weight`.
- **Group speed alignment** — the wrapper scrapes group speeds from the .prt
  file and re-interpolates them against kᵣ; this port keeps the per-mode
  values directly (equivalent). As in the Fortran, the real-path KRAKEN
  leaves VG unset (its assignment is commented out upstream), so
  `complex_solver=false` reports zeros.
- **Deterministic restarts** — KRAKENC's random restart of the secant search
  (`RANDOM_NUMBER`) is replaced by a deterministic xorshift generator seeded
  per solve, so results are reproducible run-to-run.
- **Multithreading** — inverse iteration (independent per mode) and field
  summation (independent per range column) are optionally chunked across
  threads (`threads` kwarg, default `Threads.nthreads()`). Each mode/column
  is written by exactly one task, so threaded results are bit-identical to
  serial.
- **Density subtabulation** — rho is interpolated rho-linearly in both the
  C-linear and spline SSP paths (the Fortran splines rho in `cCubic`; the
  wrapper only ever writes piecewise-constant rho per medium, for which both
  are exact).
- `attenuation.jl` is a duplicate of `src/bellhopjl/attenuation.jl` (both
  port misc/AttenMod.f90); duplicated to keep each GPL submodule
  self-contained.
- **No modefile-based initial guesses** ('I' option) and no `Solve3` profile
  continuation — single-profile scope.
- **Errors instead of dummy mode files** — "No modes for given phase speed
  interval" raises a Julia error (the Fortran writes a dummy MODFile for
  FIELD3D and calls ERROUT).
- **ForwardDiff** — not supported in this initial port (the eigenvalue
  search is not dual-safe); use the fallback `FiniteDifferences` route if
  gradients are needed.

## 4. Validation provenance

- `test/test_krakenjl.jl` cross-checks against
  `UnderwaterAcoustics.PekerisModeSolver` (mode count, kᵣ to 1e-3, mean TL
  < 1 dB) on the Pekeris scenario of `test_kraken.jl`, and at runtime against
  the Fortran `Kraken` wrapper (multilayer elastic seabed and Munk spline
  SSP: median |ΔTL| < 0.1 dB, coherent and incoherent).
- **Golden references** (hardcoded constants in `test_krakenjl.jl`) were
  generated from the Fortran Kraken/KrakenC (`AcousticsToolbox_jll`
  2025.9.6 = OALIB 2024_12_25) on 2026-07-13: Pekeris kᵣ (real to 1e-7,
  imag to 1%) and TL (0.1 dB); lossy elastic seabed kᵣ and TL likewise.
- Observed agreement with the Fortran on the scenarios above: mean |ΔTL|
  ≈ 0.003 dB (Pekeris), ≈ 0.0004 dB (elastic), ≈ 0.002 dB (Munk spline) —
  residuals dominated by the Fortran's single-precision file records.
- The real-solver (`complex_solver=false`) elastic-seabed case fails with
  "No modes for given phase speed interval" exactly as the Fortran KRAKEN
  does (it caps cHigh at the halfspace shear speed).
