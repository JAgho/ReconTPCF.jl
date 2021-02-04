
push!(LOAD_PATH,"../src/")
using Documenter, ReconTPCF, Pkg

#Pkg.add(path="C:/Users/James/.julia/dev/ReconTPCF")
makedocs(root ="./docs",
        source = "src",
        sitename="ReconTPCF.jl",
        doctest = true,
        modules = [ReconTPCF],
        pages = ["index.md", "histrecon.md", "library.md"],
        format = Documenter.HTML(prettyurls = false)
        )
