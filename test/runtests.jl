using SDPAFamily
using Test

# Make sure we get the MOI branch of Convex. This can be removed once Convex.jl proper supports MOI.
using Pkg
Pkg.add(PackageSpec(name="Convex", url="https://github.com/ericphanson/Convex.jl", rev="MathOptInterface2"))


@testset "General utilities" begin
    include("presolve.jl")
    include("variant_test.jl")
end

@testset "SDPA-GMP" begin
    global var = :gmp
    include("MOI_wrapper.jl")
    include("Convex.jl")
end

@testset "SDPA-DD" begin
    global var = :dd
    include("MOI_wrapper.jl")
    include("Convex.jl")
end

@testset "SDPA-QD" begin
    global var = :qd
    include("MOI_wrapper.jl")
    include("Convex.jl")
end

@testset "High-precision example" begin
    include("high_precision_test.jl")
end
