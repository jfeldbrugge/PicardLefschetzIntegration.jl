# PicardLefschetzIntegration

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jfeldbrugge.github.io/PicardLefschetzIntegration.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jfeldbrugge.github.io/PicardLefschetzIntegration.jl/dev/)
[![Build Status](https://github.com/jfeldbrugge/PicardLefschetzIntegration.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/jfeldbrugge/PicardLefschetzIntegration.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/jfeldbrugge/PicardLefschetzIntegration.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/jfeldbrugge/PicardLefschetzIntegration.jl)

This package numerically evaluates multidimensional oscillatory integrals of the form 
```math
I = \int e^{i S(x)}\mathrm{d}x\,,
```
for meromorphic functions $S$ using Picard Lefschetz theory. We deform the original integration domain in the complex plane onto the Lefschetz thimbles and numerically evaluate the integral.
  
## Installation
The PicardLefschetzIntegration package can be installed with the Julia package manager.
From the Julia REPL, type `]` to enter the Pkg REPL mode and run:

```julia
pkg> add PicardLefschetzIntegration
```

Or, equivalently, via the `Pkg` API:

```julia
julia> import Pkg; Pkg.add("PicardLefschetzIntegration")
```

## Usage
Please have a look at the [Tutorial page](https://jfeldbrugge.github.io/PicardLefschetzIntegration.jl/dev/) for details on how to use this package.

## Citation
When using this package for scientific research, please cite the Zenodo link [](). You can use the Bibtex entry


## Contributors
This code was written by:
* Job Feldbrugge ([job.feldbrugge@ed.ac.uk](mailto:job.feldbrugge@ed.ac.uk))
