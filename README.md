[![CI](https://github.com/org-arl/AcousticsToolbox.jl/workflows/CI/badge.svg)](https://github.com/org-arl/AcousticsToolbox.jl/actions)
[![Codecov](https://codecov.io/gh/org-arl/AcousticsToolbox.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/org-arl/AcousticsToolbox.jl)

# AcousticsToolbox

This package provides a Julia wrapper to the [OALIB](http://oalib.hlsresearch.com/AcousticsToolbox/) acoustic propagation modeling toolbox,
making it available for use with [`UnderwaterAcoustics.jl`](https://github.com/org-arl/UnderwaterAcoustics.jl).

Currently, only two of the OALIB models are supported:

- Bellhop 2D Gaussian beam tracer (almost complete support)
- Kraken 2D normal mode model (partial support)

---

## Installation

```julia
julia> # press ]
pkg> add UnderwaterAcoustics
pkg> add AcousticsToolbox
pkg> # press BACKSPACE
julia> using UnderwaterAcoustics
julia> using AcousticsToolbox
julia> models()
3-element Vector{Any}:
 PekerisRayModel
 Bellhop
 Kraken
```

## Usage

The propagation modeling API is detailed in the [UnderwaterAcoustics](https://org-arl.github.io/UnderwaterAcoustics.jl/stable/) documentation.
We assume that the reader is familiar with it. This documentation only provides guidance on specific use of `Bellhop` and `Kraken` propagation models.

## Bellhop

Additional options available with `Bellhop`:

- `nbeams` -- number of beams used for ray tracing (default: auto)
- `minangle` -- minimum beam angle in radians (default: -80°)
- `maxangle` -- maximum beam angle in radians (default: 80°)
- `gaussian` -- geometric rays if `false`, Gaussian beams if `true` (default: `false`)
- `debug` -- if `true`, intermediate Bellhop files are made available for user inspection (default: `false`)

**Example:**

```julia
using UnderwaterAcoustics
using AcousticsToolbox
using Plots

env = UnderwaterEnvironment(
  seasurface = Vacuum,
  seabed = SandyClay,
  ssp = SampledSSP(0.0:20.0:40.0, [1540.0, 1510.0, 1520.0], :smooth),
  bathymetry = SampledDepth(0.0:50.0:100.0, [40.0, 35.0, 38.0], :linear)
)
pm = Bellhop(env; gaussian=true)
tx = AcousticSource(0.0, -5.0, 1000.0)
rx = AcousticReceiverGrid2D(1.0, 0.1, 1000, -40.0, 0.2, 200)
x = transmissionloss(pm, tx, rx)
plot(env; receivers=rx, transmissionloss=x)
```

![](https://raw.githubusercontent.com/org-arl/AcousticsToolbox.jl/main/docs/images/txloss2.png)

For more information on how to use the propagation models, see [Propagation modeling toolkit](https://org-arl.github.io/UnderwaterAcoustics.jl/stable/pm_basic.html).

## Kraken

Additional options available with `Kraken`:

- `nmodes` -- maximum number of modes (default: 9999)
- `nmedia` -- number of medium (default: 1)
- `nmesh` -- number of mesh point to use initially (0=auto, default: 0)
- `clow` -- lower phase speed limit in m/s (0=auto, default: 0)
- `chigh` -- higher phase speed limit in m/s (larger values => more modes, default: 1600)
- `debug` -- if `true`, intermediate Kraken files are made available for user inspection (default: `false`)

**Example:**

```julia
using UnderwaterAcoustics
using AcousticsToolbox
using Plots

env = UnderwaterEnvironment(
  seasurface = Vacuum,
  seabed = SandyClay,
  ssp = SampledSSP(0.0:15.0:30.0, [1447.0, 1455.0, 1460.0], :smooth),
  bathymetry =  ConstantDepth(30.0)
)
pm = Kraken(env)
tx = AcousticSource(0.0, -5.0, 1000.0)
rx = AcousticReceiverGrid2D(1.0, 0.1, 5000, -30.0, 0.2, 150)
x = transmissionloss(pm, tx, rx)

plot(env; receivers=rx, transmissionloss=x, clims = (-60.0,0.0))
```

![](https://raw.githubusercontent.com/org-arl/AcousticsToolbox.jl/main/docs/images/txloss3.png)

For more information on how to use the propagation models, see [Propagation modeling toolkit](https://org-arl.github.io/UnderwaterAcoustics.jl/stable/pm_basic.html).
