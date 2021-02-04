# ReconTPCF.jl

ReconTPCF is a toolbox for the (re)construction of binary images based on the method of
 [**Torquato and Yeong**](https://doi.org/10.1103/PhysRevE.57.495). It also collects
 a set of efficient functions for computing ``S_{2}`` and ``C_{2}``, the two-point
 correlation function, and equivalent cluster correlation function. Under specialised
 conditions, functions exist to update S2 and C2 as pairs or groups of pixels change
 value.

 ## Installation

 ReconTPCF is a Julia package and hence may be installed simply using the package manager:
```@repl index
julia> ]
pkg> add https://github.com/JAgho/ReconTPCF.jl
```
This will add the package and any dependencies it has.

## Using ReconTPCF For Reconstruction

it can be run with ```using ReconTPCF``
