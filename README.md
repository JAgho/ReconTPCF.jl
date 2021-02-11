# ReconTPCF

[![Build Status](https://travis-ci.com/JAgho/ReconTPCF.jl.svg?branch=master)](https://travis-ci.com/JAgho/ReconTPCF.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/JAgho/ReconTPCF.jl?svg=true)](https://ci.appveyor.com/project/JAgho/ReconTPCF-jl)
[![Coverage](https://coveralls.io/repos/github/JAgho/ReconTPCF.jl/badge.svg?branch=master)](https://coveralls.io/github/JAgho/ReconTPCF.jl?branch=master)

Statistical Descriptor based image reconstruction of binary images. Allows for a binary image to evolve via simulated annealing into a energy minimised analogue of the original image, sharing its S2 and C2 to a major extent. Uses surface optimisation to assist the closing of shapes (limiting the range of pixels which are operated on). This library provides simple functions to compute the two-point correlation function and two point cluster correlation function in Julia. 

S2 and C2 are computed using a gridded approach, not a KD tree or ball tree. While this is detrimental to how the algorithm scales with large datasets, for smaller computations (self comparisons up to ```10^5``` points) this approach is comparable. This is achieved with fast, multithreaded computation of distance metrics (which may be arbitrary, unlike a KD tree). This software is intended to be simple, available and quick.

## Install
ReconTPCF is implemented as a full Julia package. It may be installed directly from the Julia REPL by first opening the package manager with ``]``. Then call ``Pkg add https://github.com/JAgho/ReconTPCF.jl`` and the files will be downloaded into your julia/dev folder. The toolbox contains some example executions in ``./test``, where the individual files act to load the package with ``using``. These file may simply be run using ``julia --threads 12 ./julia/dev/ReconTPCF/test/file.jl`` to use 12 cores. Alternatively, this may be run from VSCode, Atom or the REPL itself.

## Use
The test files show how the package may be called with ``using ReconTPCF``, giving access to some of the exposed functions. If internal functions are to be called, either add them to the exports in ReconTPCF.jl, or call them as a method of the package, e.g. as ```ReconTPCF.foo()```. Most functions in this toolbox are designed to work with arrays of ```CartesianIndices``` (for a basis of dimensionality N, a tuple of N Float64).

A test file may be run either from a shell by:
```path/to/.julia/dev/ReconTPCF/src/file.jl```
Or may be opened in an editor of your choosing and executed at will.

## Intended Functionality
Several modes of operation are due to be added to ReconTPCF. These will improve performance, add new ways of calculating correlation statistics, or extend the software to 3D as well as 2D data. Some possible actions are:

- Add KD tree functionality for larger datasets and profile to determine efficiency cutoff for this approach.
- Extend the approach to 3D spaces (this should be relatively simple for S2, though for C2 some clustering algorithms must be re-written).
- Add a framework for ungridded calculations of S2.
- Add k-means clustering for ungridded calculation of C2.
- Add work-order scheduling via MPI. This would send proposed pixel swaps to multiple nodes for calculation, then return the cost function of the swap.

## Testing and Compatibility
ReconTPCF is tested for VSCode (with the Julia plugin), direct REPL use, and for Juno (within Atom). It has been tested in Jupyter, but it is not determined if it is fully stable.

This package moves a lot of data, and hence task scheduling is central to its activity. This is dominated by ```Threads.@threads``` which parallelises instructions. The long nature of calculations can mean a ```@sync``` lock by another Julia program can cause the activity to be registered as a hang. Consequently, if real-time elements are included, an @async / @sync block should be provided to prevent execution while an update occurs. This is critical if e.g. ImageView.jl is to be used.
