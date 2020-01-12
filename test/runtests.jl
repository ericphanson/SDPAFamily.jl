using SDPAFamily
using Test

# Make sure we get the master branch of Convex. This can be removed once a Convex.jl release supports MOI.
using Pkg
Pkg.add(PackageSpec(name="Convex", url="https://github.com/JuliaOpt/Convex.jl", rev="master"))

const variants = (:sdpa, :sdpa_dd, :sdpa_qd, :sdpa_gmp)

@testset "SDPAFamily" begin

    @info "Starting testset `General utilities`"
    @testset "General utilities" begin
        include("presolve.jl")
        include("status_test.jl")
        include("variant_test.jl")
        include("attributes.jl")
    end

    include("Convex.jl")

    include("MOI_wrapper.jl")

    @info "Starting testset `High-precision example`"
    @testset "High-precision example" begin
        include("high_precision_test.jl")
    end
end
