using SDPAFamily
using Test

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
