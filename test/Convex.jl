using Convex
using Convex: DotMultiplyAtom
using Test
using SDPA_GMP
using MathOptInterface
using Random
const MOI = MathOptInterface
using GenericLinearAlgebra
using LinearAlgebra
import Random.shuffle
import Statistics.mean

TOL = 1e-3
eye(n) = Matrix(1.0I, n, n)

# Seed random number stream to improve test reliability
Random.seed!(2)

solvers = []

push!(solvers, () -> SDPA_GMP.Optimizer{BigFloat}(presolve = false, silent = true, variant = var))


@testset "Convex" begin
    include(joinpath("Convex", "test_utilities.jl"))
    include(joinpath("Convex", "test_const.jl"))
    include(joinpath("Convex", "test_affine.jl"))
    include(joinpath("Convex", "test_lp.jl"))
    solvers[1] = () -> SDPA_GMP.Optimizer{BigFloat}(presolve = true, silent = true, variant = var)
    include(joinpath("Convex", "test_socp.jl"))
    include(joinpath("Convex", "test_sdp.jl"))
    # include(joinpath("Convex", "test_exp.jl")
    # include(joinpath("Convex", "test_sdp_and_exp.jl")
    # include(joinpath("Convex", "test_mip.jl")
end

