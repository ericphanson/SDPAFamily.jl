using Test
using Convex
using MathOptInterface
using SparseArrays
const MOI = MathOptInterface
const MOIT = MOI.Test
const MOIU = MOI.Utilities
const MOIB = MOI.Bridges

import SDPA_GMP

@testset "Presolve: reduce subroutine" begin
    setprecision(256)
    for i in 1:10
        M = sprand(BigFloat, 10, 20000, 0.003)
        M[:, 20] .= big"0.0" # reduce does not consider the last column
        dropzeros!(M)
        j = rand(setdiff(1:10,i))
        k = rand(setdiff(1:10,i))
        λ = rand(BigFloat)
        M[i, :] = λ*M[j, :] + (1-λ)*M[k, :]
        rows = Set(rowvals(SDPA_GMP.reduce!(M)[:, 1:end-1]))
        redundant = collect(setdiff!(Set(1:10), rows))
        @test length(redundant) == 1
        @test redundant[1] ∈ (i, j, k)
    end
    for i in 1:10
        M = sprand(BigFloat, 10, 2000, 0.01)
        M[:, 20] .= big"0.0" # reduce does not consider the last column
        dropzeros!(M)
        j = rand(setdiff(1:10,i))
        k = rand(setdiff(1:10,i))
        l = rand(setdiff(1:10,i))
        λ = rand(BigFloat)
        M[i, :] = λ*M[j, :] + (1-λ)*M[k, :] + λ/2*M[l, :]
        rows = Set(rowvals(SDPA_GMP.reduce!(M)[:, 1:end-1]))
        redundant = collect(setdiff!(Set(1:10), rows))
        @test length(redundant) == 1
        @test redundant[1] ∈ (i, j, k, l)
    end
end

@testset "Presolve: redundant constraints" begin
    for n in 2:5
        opt = SDPA_GMP.Optimizer(silent = true)
        opt.no_solve = true
        A = rand(n,n) + im*rand(n,n)
        A = A + A' # now A is hermitian
        x = ComplexVariable(n,n)
        objective = sumsquares(A - x)
        c1 = x in :SDP
        problem = Problem{BigFloat}(:minimize, objective, c1)

        @test_throws BoundsError solve!(problem, opt)
        @test length(SDPA_GMP.presolve(opt)) == n^2 - n
    end
end

@testset "Presolve: inconsistent constraints" begin
    m = SDPA_GMP.Optimizer()
    m.blockdims = [3]
    m.elemdata = [(1, 1, 1, 1, big"1.0"), (1, 1, 2, 2, big"1.0"), (1, 1, 3, 3, big"1.0"),
           (2, 1, 1, 2, big"2.0"), (2, 1, 2, 1, big"2.0"),
           (3, 1, 1, 1, big"2.0"), (3, 1, 2, 2, big"2.0"), (3, 1, 3, 3, big"2.0"),
           (3, 1, 1, 2, big"3.0"), (3, 1, 2, 1, big"3.0")]
    m.b = [big"2.2", big"3.1", big"9.05"]
    @test length(SDPA_GMP.presolve(m))==1
    m.b = [big"2.2", big"3.1", big"9.05"-eps(9.05)]
    @test_throws SDPA_GMP.PresolveError SDPA_GMP.presolve(m)
end
