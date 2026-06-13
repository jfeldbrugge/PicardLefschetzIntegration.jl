```@meta
CurrentModule = PicardLefschetzIntegration
```

# PicardLefschetzIntegration

Documentation for [PicardLefschetzIntegration](https://github.com/jfeldbrugge/PicardLefschetzIntegration.jl). This package numerically evaluates multidimensional oscillatory integrals of the form 
```math
I = \int e^{i S(\bm{x})}\mathrm{d}\bm{x}\,,
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
Please have a look at the Tutorial page for details on how to use this package.

## Contributors
This code was written by:
* Job Feldbrugge ([job.feldbrugge@ed.ac.uk](mailto:job.feldbrugge@ed.ac.uk))