# PicardLefschetzIntegration

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jfeldbrugge.github.io/PicardLefschetzIntegration.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jfeldbrugge.github.io/PicardLefschetzIntegration.jl/dev/)
[![Build Status](https://github.com/jfeldbrugge/PicardLefschetzIntegration.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/jfeldbrugge/PicardLefschetzIntegration.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/jfeldbrugge/PicardLefschetzIntegration.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/jfeldbrugge/PicardLefschetzIntegration.jl)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20683184.svg)](https://doi.org/10.5281/zenodo.20683184)



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
Please have a look at the [Tutorial page](https://jfeldbrugge.github.io/PicardLefschetzIntegration.jl/stable/) for details on how to use this package.

## Citation
When using this package, please cite the paper *Oscillatory path integrals for radio astronomy* (2023) by Feldbrugge, Pen and Turok and the Zenodo link [DOI:10.5281/zenodo.20683184](https://doi.org/10.5281/zenodo.20683184). You can use the Bibtex entries

```
@ARTICLE{2023AnPhy.45169255F,
       author = {{Feldbrugge}, Job and {Pen}, Ue-Li and {Turok}, Neil},
        title = "{Oscillatory path integrals for radio astronomy}",
      journal = {Annals of Physics},
     keywords = {Interference, Lensing, Wave optics, Kirchhoff-Fresnel integral},
         year = 2023,
        month = apr,
       volume = {451},
          eid = {169255},
        pages = {169255},
          doi = {10.1016/j.aop.2023.169255},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2023AnPhy.45169255F},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}

@software{job_feldbrugge_2026_20683184,
  author       = {Job Feldbrugge},
  title        = {jfeldbrugge/PicardLefschetzIntegration.jl:
                   PicardLefschetzIntegration.jl
                  },
  month        = jun,
  year         = 2026,
  publisher    = {Zenodo},
  version      = {v1.0.0},
  doi          = {10.5281/zenodo.20683184},
  url          = {https://doi.org/10.5281/zenodo.20683184},
  swhid        = {swh:1:dir:fc800cc333cbf38115cf7296b244a966216bc821
                   ;origin=https://doi.org/10.5281/zenodo.20683183;vi
                   sit=swh:1:snp:76f139f585a46656e9c20b2ea6fb9f47b39c
                   a968;anchor=swh:1:rel:9d817d5fe15ec8731bb7ced50c6d
                   46e824306ad4;path=jfeldbrugge-
                   PicardLefschetzIntegration.jl-ec297a8
                  },
}
```

## Contributors
This code was written by:
* Job Feldbrugge ([job.feldbrugge@ed.ac.uk](mailto:job.feldbrugge@ed.ac.uk))
