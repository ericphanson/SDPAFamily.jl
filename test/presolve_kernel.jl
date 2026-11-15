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

    @testset "bounds and magnitude ranges of the variables" begin
        H = BigFloat
        none, all = typemin(Int), typemax(Int)
        # Rows `a x_j + b ≥ 0`, one variable per case:
        # free | x ≥ 0 | x ≤ 5 | x ≥ 3 | x ≤ -3 | 0 ≤ x ≤ 2^-30 | 2^-30 ≤ x ≤ 2^30 |
        # -8 ≤ x ≤ 1 | 4 ≤ x ≤ 4 | 0 ≤ x ≤ 0 | 2 ≤ x ≤ 1 (contradictory) |
        # x ≤ 6 and x ≤ 2 and x ≥ 0 | 2^-40 x ≥ 0 and 2^-40 (1 - x) ≥ 0
        t = H(2)^-40
        rows = [(2, 1, 0), (3, -1, 5), (4, 1, -3), (5, -1, -3), (6, 1, 0), (6, -1, H(2)^-30),
                (7, 1, -H(2)^-30), (7, -1, H(2)^30), (8, 1, 8), (8, -1, 1), (9, 1, -4), (9, -1, 4),
                (10, 1, 0), (10, -1, 0), (11, 1, -2), (11, -1, 1), (12, -1, 6), (12, -1, 2), (12, 1, 0),
                (13, t, 0), (13, -t, t)]
        A = sparse(1:length(rows), [r[1] for r in rows], H[r[2] for r in rows], length(rows), 13)
        b = H[r[3] for r in rows]
        sets = Any[MOI.Nonnegatives(length(rows))]
        lower, upper = SDPAFamily.variable_bounds(sets, A, b)
        @test lower == [nothing, 0, nothing, 3, nothing, 0, H(2)^-30, -8, 4, 0, 2, 0, 0]
        @test upper == [nothing, nothing, 5, nothing, -3, H(2)^-30, H(2)^30, 1, 4, 0, 1, 2, 1]
        lo, hi = SDPAFamily.magnitude_ranges(sets, A, b)
        @test lo == [none, none, none, 1, 1, none, -30, none, 2, none, none, none, none]
        @test hi == [all, all, all, all, all, -29, 31, 4, 3, all, all, 2, 1]

        # Bounds implied by rows with several variables: x1, x2 ≥ 0 and x1 + 2 x2 ≤ 8
        # give x1 ≤ 8 and x2 ≤ 4; then x3 ≥ x1 + x2 - 20 gives x3 ≥ -20, and x3 ≤ x1
        # gives x3 ≤ 8. x4 + x5 ≥ 1 gives nothing.
        A = sparse(H[1 0 0 0 0; 0 1 0 0 0; -1 -2 0 0 0; -1 -1 1 0 0; 1 0 -1 0 0; 0 0 0 1 1])
        b = H[0, 0, 8, 20, 0, -1]
        lower, upper = SDPAFamily.variable_bounds(Any[MOI.Nonnegatives(6)], A, b)
        @test lower == [0, 0, -20, nothing, nothing]
        @test upper == [8, 4, 8, nothing, nothing]

        # The diagonal entries of a PSD block are nonnegative: [x1 x2; x2 1 - x1] ⪰ 0
        # gives 0 ≤ x1 ≤ 1 and nothing for x2.
        A = sparse(H[1 0; 0 1; -1 0])
        b = H[0, 0, 1]
        lower, upper = SDPAFamily.variable_bounds(Any[MOI.PositiveSemidefiniteConeTriangle(2)], A, b)
        @test lower == [0, nothing]
        @test upper == [1, nothing]
    end

    @testset "scaling with the ranges of the variables" begin
        H = BigFloat
        scale(A, b, c) = SDPAFamily.scaled_sdpa_model(H, 256, sparse(H.(A)), H.(b), H.(c),
                                                      Any[MOI.Nonnegatives(length(b))], Int[], Int[])
        # 0 ≤ x ≤ 2^-30 next to 1 - 2^-10 x ≥ 0: equilibration and the right-hand side
        # 1 give x the unit 2^10. The range gives the unit 2^-29, and the bound row is
        # normalized again.
        s, U = big(2.0)^-10, big(2.0)^-30
        _, scaling, (As, bs, cs) = scale([1; -1; -s;;], [0, U, 1], [-1])
        @test scaling.col == [2U]
        @test As[1:2, 1] == [1, -1]
        @test bs == [0, 1 / 2, 1]
        @test abs.(cs) == [1]
        # Of the two bounds x ≤ 2^30 and x ≤ 3, the smaller one gives the range, although
        # the right-hand side of the other is the largest.
        _, scaling, _ = scale([1; -1; -1;;], [0, big(2.0)^30, 3], [-1])
        @test scaling.col == [4]
        # A variable that is at least 2^40 has at least that unit: x1 ≥ 2^40 next to
        # 2^50 - x2 ≥ 0 and x1 - x2 ≥ 0.
        _, scaling, _ = scale([1 0; 0 -1; 1 -1], [-big(2.0)^40, big(2.0)^50, 0], [1, 1])
        @test scaling.col[1] >= big(2.0)^40
    end
end
