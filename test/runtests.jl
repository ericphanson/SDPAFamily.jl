using SDPA_GMP
using Test
using Pkg

Pkg.add(PackageSpec(name="Convex", url="https://github.com/ericphanson/Convex.jl", rev="MathOptInterface"))

include("MOI_wrapper.jl")
include("Convex.jl")
include("presolve.jl")
