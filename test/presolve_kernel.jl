# Unit tests of the eliminations and of the scaling of the presolve (no solver calls).

# An integer matrix with planted row dependencies, a right-hand side `b = -A x0`, which
# is inconsistent with probability 0.2, and an objective in the row space of `A`,
# which is inconsistent with probability 0.3.
function planted_system(rng)
    m, n = rand(rng, 1:12), rand(rng, 1:12)
    A = rand(rng, -3:3, m, n) .* (rand(rng, m, n) .< rand(rng, (0.3, 0.7, 1.0)))
    for _ in 1:rand(rng, 0:3)
        i = rand(rng, 1:m)
        A[i, :] = sum(rand(rng, -2:2) * A[k, :] for k in 1:m if k != i; init = zeros(Int, n))
    end
    b = -A * rand(rng, -3:3, n)
    rand(rng) < 0.2 && (b[rand(rng, 1:m)] += 1)
    c = transpose(A) * rand(rng, -3:3, m)
    rand(rng) < 0.3 && (c[rand(rng, 1:n)] += 1)
    return A, b, c
end

@testset "Presolve: eliminations and scaling" begin

    @testset "agreement with exact arithmetic, $T" for T in (Float64, BigFloat)
        # The decisions (rank, which rows and columns are dependent, consistency) are
        # those of exact rational arithmetic; the pivots may differ.
        rng = Random.MersenneTwister(1)
        Q = Rational{BigInt}
        for _ in 1:100
            A, b, c = planted_system(rng)
            m, n = size(A)
            exact = SDPAFamily.eliminate_equalities(sparse(Q.(A)), Q.(b))
            elim = SDPAFamily.eliminate_equalities(sparse(T.(A)), T.(b))
            @test !elim.failed
            @test (elim.infeasible === nothing) == (exact.infeasible === nothing)
            @test elim.independent == exact.independent
            if elim.infeasible === nothing
                x = SDPAFamily.expand_solution(elim, n, T.(rand(rng, -3:3, length(elim.free))))
                @test norm(A * x + b, Inf) <= 1e-9 * (1 + norm(x, Inf))
            else
                @test norm(transpose(A) * elim.infeasible, Inf) <= 1e-9 * norm(elim.infeasible, Inf)
                @test dot(b, elim.infeasible) < 0
            end
            exact = SDPAFamily.find_redundant_columns(sparse(Q.(A)), Q.(c))
            result = SDPAFamily.find_redundant_columns(sparse(T.(A)), T.(c))
            @test !result.failed
            @test result.redundant == exact.redundant
            @test (result.ray === nothing) == (exact.ray === nothing)
            if result.ray !== nothing
                @test norm(A * result.ray, Inf) <= 1e-9 * norm(result.ray, Inf)
                @test dot(c, result.ray) < 0
            end
        end
    end

    @testset "decisions do not depend on the scaling, $T" for T in (Float64, BigFloat)
        # Rows and columns scaled by powers of two (exact): the same rows and columns
        # are dependent.
        rng = Random.MersenneTwister(2)
        kmax = T == Float64 ? 100 : 5000
        for _ in 1:50
            A, b, c = planted_system(rng)
            m, n = size(A)
            r = T(2) .^ rand(rng, -kmax:kmax, m)
            s = T(2) .^ rand(rng, -kmax:kmax, n)
            As = sparse(Diagonal(r) * T.(A) * Diagonal(s))
            plain = SDPAFamily.eliminate_equalities(sparse(T.(A)), T.(b))
            scaled = SDPAFamily.eliminate_equalities(As, r .* T.(b))
            @test !scaled.failed
            @test scaled.independent == plain.independent
            @test (scaled.infeasible === nothing) == (plain.infeasible === nothing)
            plain = SDPAFamily.find_redundant_columns(sparse(T.(A)), T.(c))
            scaled = SDPAFamily.find_redundant_columns(As, s .* T.(c))
            @test !scaled.failed
            @test scaled.redundant == plain.redundant
            @test (scaled.ray === nothing) == (plain.ray === nothing)
        end
    end

    @testset "dense system with a pivot order driven by sparsity, $T" for T in (Float64, BigFloat)
        # Random pivot costs make the elimination prefer pivots down to 1/10 of the
        # largest entry. The elimination must not lose accuracy with the size.
        rng = Random.MersenneTwister(3)
        m, n = 150, 170
        A = sparse(T.(randn(rng, m, n)))
        b = -A * T.(randn(rng, n))
        elim = SDPAFamily.eliminate_equalities(A, b; pivot_cost = rand(rng, 1:5, n))
        @test !elim.failed
        @test length(elim.pivots) == m
        x = SDPAFamily.expand_solution(elim, n, T.(randn(rng, n - m)))
        @test maximum(abs.(A * x + b) ./ (abs.(A) * abs.(x) + abs.(b))) <= 100 * eps(T)
        g = transpose(A) * T.(randn(rng, m))
        y = SDPAFamily.equality_duals(elim, m, g)
        @test norm(transpose(A) * y - g, Inf) <= 1000 * eps(T) * norm(transpose(abs.(A)) * abs.(y), Inf)
    end

    @testset "rounding residue and small exact data, $T" for T in (Float64, BigFloat)
        # Row 2 is 3 * row 1 plus the unit vector of variable 3. The elimination
        # divides by 3; the residue of the cancellation must not become an entry.
        AE = sparse(T[1 1 0 1; 3 3 1 3] ./ 3)
        elim = SDPAFamily.eliminate_equalities(AE, T[0, 0])
        @test length(elim.pivots) == 2
        @test sum(length(row.idx) for row in elim.rows) == 2
        # x == y and L x - L y + s w == 1: the coefficient s is exact data and stays,
        # although it is far below L (but not below the rounding residue of L).
        L, s = T == Float64 ? (1e8, 1e-8) : (big"1e40", big"1e-40")
        AE = sparse(T[1 -1 0; L -L s])
        elim = SDPAFamily.eliminate_equalities(AE, T[0, -1])
        @test elim.infeasible === nothing
        @test length(elim.pivots) == 2
        x = SDPAFamily.expand_solution(elim, 3, T[5])
        @test x[3] ≈ 1 / s
        @test x[1] == x[2]
    end

    @testset "scaling factors" begin
        rng = Random.MersenneTwister(4)
        # a matrix with entries between 1 and 4, rows and columns scaled by 2^±300
        N = sprand(rng, 30, 20, 0.3) .* 3
        N.nzval .+= 1
        r, s = rand(rng, -300:300, 30), rand(rng, -300:300, 20)
        M = Diagonal(big(2.0) .^ r) * BigFloat.(N) * Diagonal(big(2.0) .^ s)
        ge, se = SDPAFamily.equilibrate(M, collect(1:30), 30)
        scaled = Diagonal(big(2.0) .^ ge) * M * Diagonal(big(2.0) .^ se)
        for i in 1:30
            any(!iszero, N[i, :]) && @test 1 / 2 <= maximum(abs, scaled[i, :]) < 2
            any(!iszero, N[i, :]) || @test ge[i] == 0
        end
        for j in 1:20
            any(!iszero, N[:, j]) && @test 1 / 2 <= maximum(abs, scaled[:, j]) < 2
        end
        # rows of one group share a factor
        ge, se = SDPAFamily.equilibrate(M, repeat(1:10, inner = 3), 10)
        @test length(ge) == 10
        # Two parts that share no variable, at different scales:
        # -x1 + 2^20 >= 0, 2^-10 x2 + 2^31 >= 0, min -2^-18 x1 + 2^-39 x2.
        # Each part gets its own factor for the right-hand side.
        A = sparse(BigFloat[-1 0; 0 big(2.0)^-10])
        b = BigFloat[big(2.0)^20, big(2.0)^31]
        c = BigFloat[-big(2.0)^-18, big(2.0)^-39]
        sets = Any[MOI.Nonnegatives(2)]
        _, scaling, (As, bs, cs) = SDPAFamily.scaled_sdpa_model(BigFloat, 256, A, b, c, sets, Int[], Int[])
        @test abs.(Matrix(As)) == [1 0; 0 1]
        @test bs == [1, 1]
        @test abs.(cs) == [1, 1]
        @test scaling.row .* (A * (scaling.col .* [1, 1])) == As * [1, 1]
    end
end
