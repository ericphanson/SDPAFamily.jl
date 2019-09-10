using SDPAFamily
using Test

# Make sure we get the MOI branch of Convex. This can be removed once Convex.jl proper supports MOI.
using Pkg
Pkg.add(PackageSpec(name="Convex", url="https://github.com/ericphanson/Convex.jl", rev="MathOptInterface"))

@testset "SDPAFamily" begin

    @testset "General utilities" begin
        include("presolve.jl")
        include("variant_test.jl")
    end

    @testset "SDPA-plain" begin
        global var = :sdpa
        include("MOI_wrapper.jl")
        include("Convex.jl")
    end

    @testset "SDPA-DD" begin
        global var = :sdpa_dd
        include("MOI_wrapper.jl")
        include("Convex.jl")
    end

    @testset "SDPA-QD" begin
        global var = :sdpa_qd
        include("MOI_wrapper.jl")
        include("Convex.jl")
    end

    @testset "SDPA-GMP" begin
        global var = :sdpa_gmp
        include("MOI_wrapper.jl")
        include("Convex.jl")
    end

    @testset "High-precision example" begin
        include("high_precision_test.jl")
    end
end
