# BellhopJL porting notes — provenance and deviations

This directory contains a native Julia port of the 2D BELLHOP Gaussian beam /
ray tracer. This file documents the port's lineage, every intentional deviation
from the Fortran sources, and the provenance of the validation data.
(Everything in `src/bellhopjl/` is GPL-3.0-or-later — see `LICENSE-GPL-3.0` at
the repository root; the rest of AcousticsToolbox.jl is MIT.)

## 1. Lineage

Ported from **BELLHOP 2D as released in the OALIB Acoustics Toolbox 2022_4**,
via the **A-New-BellHope mirror** (https://github.com/A-New-BellHope/bellhop,
© 2021–2023 The Regents of the University of California; BELLHOP © 1983–2022
Michael B. Porter), *including* its documented bug fixes. Docstrings and code
comments cite the originating Fortran file/subroutine throughout (e.g. "ports
`InfluenceGeoHatCart`, influence.f90") so the code remains diffable against
upstream.

The following **A-New-BellHope deviations from stock OALIB BELLHOP are
deliberately retained** (stock OALIB, including the 2024_12_25 release, does
not have them). They fix the reproducibility problem documented in the
A-New-BellHope README — in stock BELLHOP the reduced integration step lands
"randomly" just before or just after every boundary/interface crossing, so
*every step of every ray is an edge case* whose resolution can differ between
compilers, runs and implementations:

- **Exact boundary landing (`StepToBdry2D`, Step.f90)** — the blended step is
  placed *exactly* on the nearest SSP interface / boundary / box limit, with
  snap-to-depth for flat boundaries, and reflections are signalled by explicit
  `topRefl`/`botRefl` flags rather than a distance-sign test. Retained because
  it makes results deterministic and well-defined; reverting it would
  reintroduce the documented edge-case chaos.
- **Tangent-direction-aware segment selection** (`GetTopSeg`/`GetBotSeg`,
  bdryMod.f90; `UpdateDepthSegmentT`, sspMod.f90) — after landing exactly on a
  node, the segment is chosen so that a small step along the ray stays inside
  it. Required for the exact-landing scheme to be consistent.
- **Small-step scheme** — `INFINITESIMAL_STEP_SIZE = 1e-6·deltas` (stock 2022
  used a 1e-4 threshold with a 1e-5 makeup step; 2024 uses 1e-4/1e-4). The
  smaller constant belongs with exact landing, which makes tiny cleanup steps
  rare.
- **Gaussian shallow-angle threshold 0.50001** (`InfluenceGeoGaussianCart`,
  influence.f90) — moves a behavioural discontinuity off the round number so
  rays launched at exactly 60° do not straddle it.

## 2. Subsequently adopted OALIB 2024_12_25 corrections

After the initial port, the 2D-relevant algorithmic changes of the official
OALIB release 2024_12_25 (relative to the 2022_4 base) were reviewed and the
following adopted:

- **influence.f90 : `InfluenceGeoHatCart`, `InfluenceGeoGaussianCart`** — the
  caustic phase shift at a receiver whose interpolated `q` straddles a caustic
  is now *additive* (`phaseInt = phaseInt + π/2`); the 2022 code discarded the
  accumulated ray phase (`phaseInt = phase + π/2`), a bug A-New-BellHope had
  flagged but not fixed in 2D.
- **influence.f90 : `InfluenceGeoGaussianCart`** — the ray phase is read from
  the *previous* ray point (`ray2D(iS-1)%Phase`), consistent with the hat
  influence; 2022 read the current point.
- **Step.f90 : `Step2D`** — the SSP-interface jump correction on `p` is applied
  only for SSP types with a discontinuous first derivative ('N'/'C'); C¹
  profiles (spline 'S', PCHIP 'P', analytic 'A') get `gradcjump = 0`
  (implemented as the `has_gradcjump` trait in `ssp.jl`).
- **bellhopMod.f90** — `MaxN` (per-ray step cap) raised 100 000 → 1 000 000.

Changes reviewed and **not** adopted: the 1e-4 small-step makeup constant (tied
to the stock no-snapping stepping, see §1), receiver-position/env-file I/O
precision changes (no file I/O here), Cerveny-beam changes (out of scope), and
cosmetics. PCHIP (`misc/pchipMod.f90`), the spline coefficients
(`misc/splinec.f90`), `Reflect2D` and the `AddArr` arrival merging are
unchanged between 2022_4 and 2024_12_25.

## 3. Other intentional deviations from the Fortran

- **No env/ray/arr/shd file I/O** — environments come directly from
  UnderwaterAcoustics.jl; the adapter (`api.jl`) mirrors the unit conversions
  of this package's Fortran env-file writer verbatim (density ratio in g/cm³,
  attenuation in dB/λ, Francois-Garrison volume attenuation with
  z̄ = waterdepth/2, angle sign flips, 1.01× ray box, SSP node layout).
- **Scope** — only the features reachable through the UnderwaterAcoustics.jl
  API are implemented: geometric hat/Gaussian beams in Cartesian coordinates,
  coherent/incoherent/semicoherent fields, arrivals/eigenrays. Cerveny and SGB
  beams, ray-centered influences, 3D, and range-dependent SSPs ('Q'/'H') are
  out of scope (the wrapper never selects them). Source beam patterns are not
  implemented (the wrapper never writes SBP files).
- **Semicoherent runs** are implemented via the Lloyd-mirror source amplitude
  `√2·|sin(ω·zₛ·sin α/c)|` at launch (BellhopCore), matching the Fortran.
- **PCHIP SSP not ported** — the wrapper maps interpolations only to
  'C'/'S'/'Q', so 'P' is unreachable. C-linear, N²-linear and cubic-spline
  (de Boor not-a-knot, `misc/splinec.f90` translated literally) are ported.
- **Curvilinear ('C') boundary interpolation not implemented** — boundaries are
  piecewise linear with the node-extension scheme of
  `ComputeBdryTangentNormal`. The wrapper emits 'C' only for non-`Linear()`
  `SampledFieldX` bathymetry; flagged for future work.
- **Per-step receiver bracketing** — the monotone receiver-pointer walk of
  `InfluenceGeoHatCart`/`InfluenceGeoGaussianCart` is ported faithfully; the
  Fortran `MINLOC` start-index edge case (no receiver to the right of the
  source ⇒ index 0, undefined behaviour) is guarded by clamping.
- **Precision** — double precision throughout (the Fortran accumulates the TL
  field in single-precision `COMPLEX*4` and quantizes .shd/.arr output).
- **AD requirements** — all numeric structs are parametric so
  `ForwardDiff.Dual`s flow through any single parameter class (SSP nodes,
  bathymetry, seabed ρ/c/α, source/receiver coordinates, frequency);
  `Env2D` permits heterogeneous element types across components; no `Float64`
  hard-casts or `@fastmath` on dual-reachable paths.
- The SSP interface jump is skipped for exactly horizontal rays (`t₂ = 0`
  guard; the Fortran would divide by zero).
- `AddArr` merges against the *last* arrival only, with `PhaseTol = 0.05` —
  Fortran behaviour retained for parity (flagged as a bug by A-New-BellHope).
- **Ray-history buffer reuse** — `TraceRay2D`'s ray storage is a fresh
  `Vector{RayPt}` per beam in the naive port; `trace_ray!` reuses one buffer
  across the beam loop (the Fortran uses a single static `ray2D` array, so
  this is actually *closer* to the original).
- **Multithreaded beam loop** — the beam loop of `BellhopCore` is optionally
  split into contiguous chunks across threads (`threads` kwarg of `BellhopJL`,
  default `Threads.nthreads()`). Each chunk has its own field accumulator /
  arrival sink; per-chunk fields are summed in chunk order (deterministic for
  a fixed thread count, differs from serial only by floating-point
  reassociation), and per-chunk arrival lists are replayed through `add_arr!`
  in beam order to preserve the serial merge semantics. `threads=1` runs the
  identical serial path.
- **Eigenray retrace memoization** — `arrivals` traces each *unique* take-off
  angle once when extracting ray paths (merged arrivals share angles); the
  Fortran writes ray files in a separate run and has no equivalent step.
- `@inbounds` on the influence receiver loop and SSP segment scan (bounds
  are guaranteed by the monotone pointer walk / segment clamping).

## 4. Validation provenance

- Analytic unit tests: straight rays and exact travel times (isovelocity),
  circular-arc rays in a linear gradient, reflection phase/amplitude, Rayleigh
  coefficient cross-checked against `UnderwaterAcoustics.reflection_coef`.
- **Golden references** (hardcoded constants in `test/test_bellhopjl.jl`) were
  generated with the Fortran BELLHOP shipped by **AcousticsToolbox_jll
  v2025.9.6**, which is built from the official OALIB `at_2024_12_25.zip`
  (SHA-pinned in Yggdrasil), driven through this package's `Bellhop` wrapper on
  Pekeris (coherent + incoherent TL, arrivals) and Munk spline-SSP (Gaussian
  beams, incoherent TL) scenarios. A direct in-process `BellhopJL`-vs-`Bellhop`
  consistency test asserts median |ΔTL| < 0.1 dB on a Pekeris grid.
- **Known residual**: on Munk-type spline SSPs the *coherent* TL of this port
  differs from the OALIB binary by ~1–2 dB median (amplitude/incoherent
  agreement is ~0.02 dB median). This is explained by the retained
  A-New-BellHope stepping (§1): stock OALIB lands randomly before/after each
  of the many SSP-interface crossings, accumulating O(1e-4–1e-3 s) delay
  differences on multi-bounce paths — the same order as the OALIB binary's own
  sensitivity to a ±1 change in beam count. Amplitude-level agreement is
  therefore tested tightly and coherent spline-SSP agreement statistically.
