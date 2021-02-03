# ReconTPCF

[![Build Status](https://travis-ci.com/JAgho/ReconTPCF.jl.svg?branch=master)](https://travis-ci.com/JAgho/ReconTPCF.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/JAgho/ReconTPCF.jl?svg=true)](https://ci.appveyor.com/project/JAgho/ReconTPCF-jl)
[![Coverage](https://coveralls.io/repos/github/JAgho/ReconTPCF.jl/badge.svg?branch=master)](https://coveralls.io/github/JAgho/ReconTPCF.jl?branch=master)

Statistical Descriptor based image reconstruction of binary images. Allows for a binary image to evolve via simulated annealing into a energy minimised analogue of the original image, sharing its S2 and C2 to a major extent. Uses surface optimisation to assist the closing of shapes (limiting the range of pixels which are operated on)

## Install
Simply add this as a package. Dependencies should be automatically fetched. Call ``Pkg add https://github.com/JAgho/ReconTPCF.jl`` and the files will be downloaded into your julia/dev folder.

## Use
The test files show how the package may be called with ``using ReconTPCF``, giving access to some of the exposed functions. If internal functions are to be called, either add them to the exports in ReconTPCF.jl, or call them as a method of the package, e.g. as ReconTPCF.foo()
