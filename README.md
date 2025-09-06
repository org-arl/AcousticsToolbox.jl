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

For information on how to use the models, see [documentation](https://org-arl.github.io/UnderwaterAcoustics.jl/).

## Contributing

Contributions in the form of bug reports, feature requests, ideas/suggestions, bug fixes, code enhancements, and documentation updates are most welcome. Please read [contribution guidelines](https://github.com/org-arl/UnderwaterAcoustics.jl/blob/master/CONTRIBUTING.md) if you wish to start contributing.

The scopes active in this repository are:
- **bellhop**: Bellhop
- **kraken**: Kraken
- **orca**: Orca
