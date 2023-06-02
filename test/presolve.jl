using Test
using Convex
import MathOptInterface as MOI
using SparseArrays
const MOIT = MOI.Test
const MOIB = MOI.Bridges

import SDPAFamily
using Random

Random.seed!(5)
setprecision(256)

# In the following, we use this of the random API:
# https://docs.julialang.org/en/v1/stdlib/Random/index.html#Generating-values-from-a-collection-1
struct RandomRationalType end
const RandomRational = RandomRationalType()

"""
Random.rand(rng::AbstractRNG, ::Random.SamplerTrivial{RandomRationalType})::Rational{BigInt}

For the purposes of these tests, we define `rand(RandomRational)` to give a random rational number in
the interval [1/10_000, 1 - 1/10_000] representable with a denominator of size at most 10_000, with the type `Rational{BigInt}`.
"""
function Random.rand(rng::AbstractRNG, ::Random.SamplerTrivial{RandomRationalType})
    D = big(rand(2:10_000))
    N = big(rand(1:D-1))
    return N // D
end

@testset "Presolve" begin

    @testset "reduce subroutine" begin
        for i in 1:10
            # make sprand generate its nonzeros via our `RandomRational`'s
            M = 1000*sprand(10, 20000, 0.003, n -> [ rand(RandomRational) for _ = 1:n ], Rational{BigInt})
            M[:, 20] .= big"0.0" # reduce does not consider the last column
            dropzeros!(M)
            j = rand(setdiff(1:10,i))
            k = rand(setdiff(1:10,i))
            λ = rand(RandomRational)
            M[i, :] = λ*M[j, :] + (1-λ)*M[k, :]
            rows = Set(rowvals(SDPAFamily.reduce!(M)[:, 1:end-1]))
            redundant = collect(setdiff!(Set(1:10), rows))
            @test length(redundant) == 1
            @test redundant[1] ∈ (i, j, k)
        end
        for i in 1:10
            M = 1000*sprand(10, 2000, 0.01, n -> [ rand(RandomRational) for _ = 1:n ], Rational{BigInt})
            M[:, 20] .= big"0.0" # reduce does not consider the last column
            dropzeros!(M)
            j = rand(setdiff(1:10,i))
            k = rand(setdiff(1:10,i))
            l = rand(setdiff(1:10,i))
            λ = rand(RandomRational)
            M[i, :] = λ*M[j, :] + (1-λ)*M[k, :] + λ/2*M[l, :]
            rows = Set(rowvals(SDPAFamily.reduce!(M)[:, 1:end-1]))
            redundant = collect(setdiff!(Set(1:10), rows))
            @test length(redundant) == 1
            @test redundant[1] ∈ (i, j, k, l)
        end
        for i in 1:20
            M = sprand(20, 100, 0.5)
            M[i, :] .= eps(norm(M, Inf))
            M[:, 1] .= eps(norm(M, Inf))
            M = SDPAFamily.reduce!(M)
            I, J, V = findnz(M)
            @test i ∉ I
        end

    end

    @testset "redundant constraints" begin
        for n in 2:5
            opt = SDPAFamily.Optimizer{Float64}(;silent = true, variant = :sdpa, presolve = true)
            opt.no_solve = true
            A = rand(Float64, n,n) + im*rand(Float64, n,n)
            A = A + A' # now A is hermitian
            x = ComplexVariable(n,n)
            p = minimize(sumsquares(A - x), [x in :SDP])

            @test_throws BoundsError solve!(p, opt)
            @test length(SDPAFamily.presolve(opt)) == n^2 - n
        end
    end

    @testset "inconsistent constraints" begin
        m = SDPAFamily.Optimizer(verbose = SDPAFamily.WARN)
        m.blockdims = [3]
        m.elemdata = [(1, 1, 1, 1, big"1.0"), (1, 1, 2, 2, big"1.0"), (1, 1, 3, 3, big"1.0"),
            (2, 1, 1, 2, big"2.0"), (2, 1, 2, 1, big"2.0"),
            (3, 1, 1, 1, big"2.0"), (3, 1, 2, 2, big"2.0"), (3, 1, 3, 3, big"2.0"),
            (3, 1, 1, 2, big"3.0"), (3, 1, 2, 1, big"3.0")]
        m.b = [big"2.2", big"3.1", big"9.05"]
        @test length(SDPAFamily.presolve(m))==1
        m.b = [big"2.2", big"3.1", big"9.05"-eps(9.05)]
        @test_logs (:warn, "Inconsistency at constraint index 1. Problem is dual infeasible.") SDPAFamily.presolve(m)
    end

end
