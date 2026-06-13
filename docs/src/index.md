```@meta
CurrentModule = PicardLefschetzIntegration
```

# PicardLefschetzIntegration

Documentation for [PicardLefschetzIntegration](https://github.com/jfeldbrugge/PicardLefschetzIntegration.jl). This package is particle mesh code in Julia that is inspired by Johan Hidding's [nbody2d](https://zenodo.org/records/4158731) Python code (see [jhidding.github.io/nbody2d](https://jhidding.github.io/nbody2d/) for more details). It serves as a quick environment to learn and experiment with cosmological $N$-body simulations in the two-dimensional setting. This enables one to quickly develop intuition. 

Using the package:

* We create initial conditions by sampling both unconstrained and constrained Gaussian random fields.
* We evolve the initial perturbations to the current epoch.
* We estimate the density, velocity and number of stream fields using the Phase-Space Delaunay Tessellation Field Estimator.
  
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
Please have a look at the Tutorial page for details on how to use this package.

## Contributors
This code was written by:
* Job Feldbrugge ([job.feldbrugge@ed.ac.uk](mailto:job.feldbrugge@ed.ac.uk))