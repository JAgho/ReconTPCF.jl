
push!(LOAD_PATH,"../src/")
using Documenter, ReconTPCF, Pkg

#Pkg.add(path="C:/Users/James/.julia/dev/ReconTPCF")
makedocs(root ="./docs",
        source = "src",
        sitename="My Documentation",
        doctest = true,
        modules = [ReconTPCF],
        pages = ["index.md"]
        )
