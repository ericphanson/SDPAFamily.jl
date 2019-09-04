using SDPA_GMP
using Test

# Make sure we get the MOI branch of Convex. This can be removed once Convex.jl proper supports MOI.
using Pkg
Pkg.add(PackageSpec(name="Convex", url="https://github.com/ericphanson/Convex.jl", rev="MathOptInterface"))

include("presolve.jl")
include("Convex.jl")
include("MOI_wrapper.jl")
