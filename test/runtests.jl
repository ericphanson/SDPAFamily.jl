using Convex, MathOptInterface, SDPAFamily
using Test

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
