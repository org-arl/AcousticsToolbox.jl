[![doc-stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://org-arl.github.io/UnderwaterAcoustics.jl/bellhop.html)
[![CI](https://github.com/org-arl/AcousticsToolbox.jl/workflows/CI/badge.svg)](https://github.com/org-arl/AcousticsToolbox.jl/actions)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![Codecov](https://codecov.io/gh/org-arl/AcousticsToolbox.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/org-arl/AcousticsToolbox.jl)
[![ColPrac](https://img.shields.io/badge/ColPrac-contributing-blueviolet)](https://github.com/org-arl/UnderwaterAcoustics.jl/blob/master/CONTRIBUTING.md)

# AcousticsToolbox

This package provides a Julia wrapper to the [OALIB](http://oalib.hlsresearch.com/AcousticsToolbox/) acoustic propagation modeling toolbox
(and other related tools), making it available for use with [`UnderwaterAcoustics.jl`](https://github.com/org-arl/UnderwaterAcoustics.jl).

Currently, the following models are supported:

- Bellhop 2D Gaussian beam tracer
- Kraken 2D normal mode model
- Orca 2D normal mode model
- BellhopJL — a native Julia port of the 2D Bellhop Gaussian beam tracer
  (no file I/O, ForwardDiff-differentiable; see
  [src/bellhopjl/PORTING_NOTES.md](src/bellhopjl/PORTING_NOTES.md) for the
  port's provenance and deviations from the Fortran)
- KrakenJL — a native Julia port of the Kraken/KrakenC 2D normal mode models
  (no file I/O, multithreaded; see
  [src/krakenjl/PORTING_NOTES.md](src/krakenjl/PORTING_NOTES.md) for the
  port's provenance and deviations from the Fortran)

For information on how to use the models, see [documentation](https://org-arl.github.io/UnderwaterAcoustics.jl/).

## Licensing

This package is dual-licensed on a per-directory basis:

- Everything **except** `src/bellhopjl/` and `src/krakenjl/` is licensed under
  the [MIT license](LICENSE).
- `src/bellhopjl/` (the `BellhopJL` solver) is a native Julia port of Bellhop — a
  derivative work of the GPL-licensed Bellhop (© 1983–2022 Michael B. Porter;
  [A-New-BellHope](https://github.com/A-New-BellHope/bellhop) changes © 2021–2023
  The Regents of the University of California) — and is licensed under
  [GPL-3.0-or-later](LICENSE-GPL-3.0).
- `src/krakenjl/` (the `KrakenJL` solver) is a native Julia port of
  Kraken/KrakenC — a derivative work of the GPL-licensed KRAKEN (© Michael B.
  Porter, OALIB Acoustics Toolbox 2024_12_25) — and is likewise licensed under
  [GPL-3.0-or-later](LICENSE-GPL-3.0).

Each file in those directories carries an `SPDX-License-Identifier:
GPL-3.0-or-later` header. The ports' lineage and every intentional deviation
from the Fortran are documented in
[src/bellhopjl/PORTING_NOTES.md](src/bellhopjl/PORTING_NOTES.md) and
[src/krakenjl/PORTING_NOTES.md](src/krakenjl/PORTING_NOTES.md).

Note that this package already downloads and executes the GPL-licensed OALIB
Fortran binaries (via `AcousticsToolbox_jll`) for the `Bellhop`/`Kraken` models, so
GPL software is involved either way; with the `BellhopJL` and `KrakenJL` components
included in the source tree, **distribution of the combined package is subject to
the terms of the GPL-3.0**. If you need MIT-only terms, strip `src/bellhopjl/` and
`src/krakenjl/` (the rest of the package does not depend on them).

## Contributing

Contributions in the form of bug reports, feature requests, ideas/suggestions, bug fixes, code enhancements, and documentation updates are most welcome. Please read [contribution guidelines](https://github.com/org-arl/UnderwaterAcoustics.jl/blob/master/CONTRIBUTING.md) if you wish to start contributing.

The scopes active in this repository are:
- **bellhop**: Bellhop
- **bellhopjl**: BellhopJL
- **kraken**: Kraken
- **krakenjl**: KrakenJL
- **orca**: Orca
