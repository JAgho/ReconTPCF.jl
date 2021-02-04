# ReconTPCF.jl

ReconTPCF is a toolbox for the (re)construction of binary images based on the method of
 [**Torquato and Yeong**](https://doi.org/10.1103/PhysRevE.57.495). It also collects
 a set of efficient functions for computing ``S_{2}`` and ``C_{2}``, the two-point
 correlation function, and equivalent cluster correlation function. Under specialised
 conditions, functions exist to update S2 and C2 as pairs or groups of pixels change
 value.

 ## Installation

 ReconTPCF is a Julia package and hence may be installed simply using the package manager:
```
    julia> ]
    pkg> add https://github.com/JAgho/ReconTPCF.jl
```
This will add the package and any dependencies it has.

## Using ReconTPCF For Reconstruction

Calling ```using ReconTPCF``` will give access to its
exported functions. This includes ```get_C2_S2(fname)``` and ```histrecon()``` which are an initialiser and the main reconstruction loop respectively. Histrecon is presently adjusted by modifying its function definition. A typical use of it would look like:
```
    dims, C2, S2, philen = get_C2_S2(fname)
    guess, S2n, C2n, S2_BN1, C2_BN1, SN1 = histrecon((200, 200)), C2, S2, 12000)
```
Further detail regarding the algorithm is given in its own page in the sidebar


## Using ReconTPCF for Computing The Two-Point Correlation Function

ReconTPCF has fast algorithms for computing S2 for large and small collections in several ways. For long lists of points, a multithreaded implementation is provided; ```blas_stat5```. For shorter lists a single thread function is given; ```blas_stat_st2```. These functions consider all unique pairs in the list and compute the L2 norm for these. These are histogrammed according to the secondary arguments of these functions.
