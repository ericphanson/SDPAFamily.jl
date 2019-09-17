using SDPAFamily
using Test

# Make sure we get the MOI branch of Convex. This can be removed once Convex.jl proper supports MOI.
using Pkg
Pkg.add(PackageSpec(name="Convex", url="https://github.com/ericphanson/Convex.jl", rev="MathOptInterface"))

const variants = (:sdpa, :sdpa_dd, :sdpa_qd, :sdpa_gmp)

@testset "SDPAFamily" begin

    @testset "General utilities" begin
        include("presolve.jl")
        include("variant_test.jl")
    end

    include("Convex.jl")

    include("MOI_wrapper.jl")

    @testset "High-precision example" begin
        include("high_precision_test.jl")
    end
end
