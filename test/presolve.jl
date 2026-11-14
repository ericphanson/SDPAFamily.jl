using Test
using Convex
import MathOptInterface as MOI
using SparseArrays
using LinearAlgebra

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

random_rational_matrix(m, n, p) =
    1000 * sprand(m, n, p, k -> [rand(RandomRational) for _ = 1:k], Rational{BigInt})

# A model with the optimizer behind a cache and bridges, as JuMP and Convex use it.
function bridged_model(::Type{T}; kwargs...) where {T}
    optimizer = SDPAFamily.Optimizer{T}(; silent = true, presolve = true, kwargs...)
    cache = MOI.Utilities.UniversalFallback(MOI.Utilities.Model{T}())
    cached = MOI.Utilities.CachingOptimizer(cache, optimizer)
    return MOI.Bridges.full_bridge_optimizer(cached, T), optimizer
end

# The accuracy to expect from the solution returned by `variant` with number type `T`:
# the default `epsilonStar` and `epsilonDash` are 1e-7 for `:sdpa` and for `Float64`
# with `:sdpa_gmp`, 1e-15 for `:sdpa_dd`, 1e-20 for `:sdpa_qd` and 1e-30 for `:sdpa_gmp`;
# with `Float64`, the solution is also read back in `Float64`.
function soltol(var, T)
    var == :sdpa && return 1e-6
    T == Float64 && return 1e-6
    var == :sdpa_dd && return 1e-12
    var == :sdpa_qd && return 1e-16
    return 1e-25
end

vaf(T, terms, constant) = MOI.Utilities.vectorize([sum(T(a) * x for (a, x) in terms) + T(constant)])

@testset "Presolve" begin

    @testset "find_redundant_columns, exact arithmetic" begin
        for i in 1:10
            A = random_rational_matrix(2000, 10, 0.01)
            j = rand(setdiff(1:10, i))
            k = rand(setdiff(1:10, i, j))
            λ = rand(RandomRational)
            A[:, i] = λ * A[:, j] + (1 - λ) * A[:, k]
            dropzeros!(A)
            # consistent objective: `c` is in the row space of `A`
            c = transpose(A) * [rand(RandomRational) for _ in 1:2000]
            result = SDPAFamily.find_redundant_columns(A, c)
            @test result.ray === nothing
            @test length(result.redundant) == 1
            @test result.redundant[1] ∈ (i, j, k)
            # inconsistent objective: the result is a ray `d` with `A d = 0`, `c'd < 0`
            c[i] += 1
            result = SDPAFamily.find_redundant_columns(A, c)
            @test isempty(result.redundant)
            @test iszero(A * result.ray)
            @test dot(c, result.ray) < 0
        end
    end

    @testset "find_redundant_columns, $T" for T in (Float64, BigFloat)
        A = sparse(T.(rand(-5:5, 40, 6)))
        A[:, 6] = A[:, 1] / 3 + A[:, 2] * T(7) / 10 - A[:, 5]
        c = transpose(A) * T.(rand(40))
        result = SDPAFamily.find_redundant_columns(A, c)
        @test result.redundant == [6]
        @test result.ray === nothing
        # a perturbation that is large compared to the tolerance makes the columns independent
        A[1, 6] += T(1) / 1000
        result = SDPAFamily.find_redundant_columns(A, c)
        @test isempty(result.redundant)
    end

    @testset "eliminate_equalities, $T" for T in (Rational{BigInt}, Float64, BigFloat)
        AE = T.(rand(-3:3, 4, 7))
        AE[4, :] = AE[1, :] + 2 * AE[2, :]
        AE = sparse(AE)
        x0 = T.(rand(-3:3, 7))
        bE = -AE * x0
        elim = SDPAFamily.eliminate_equalities(AE, bE)
        @test elim.infeasible === nothing
        @test length(elim.pivots) == 3
        @test length(elim.free) == 4
        # every choice of the free variables gives a solution
        x = SDPAFamily.expand_solution(elim, 7, T.(rand(-3:3, 4)))
        @test norm(AE * x + bE, Inf) <= 1e-10
        # duals: `AE' y = g` for `g` in the row space of `AE`
        g = transpose(AE) * T.(rand(-3:3, 4))
        y = SDPAFamily.equality_duals(elim, 4, g)
        @test norm(transpose(AE) * y - g, Inf) <= 1e-10
        @test iszero(y[4]) # the dependent row has dual zero
        # inconsistent right-hand side
        bE[4] += 1
        elim = SDPAFamily.eliminate_equalities(AE, bE)
        y = elim.infeasible
        @test norm(transpose(AE) * y, Inf) <= 1e-10
        @test dot(bE, y) < 0
    end

    @testset "MOI with variant $var and type $T" for var in variants, T in (Float64, BigFloat)
        tol = soltol(var, T)
        @testset "linearly dependent variables" begin
            # The columns of `x` and `y` are equal; without presolve, SDPA-DD and SDPA-QD
            # stop at the iteration limit and SDPA-GMP returns a point whose objective
            # value differs from the reported optimal value.
            model, optimizer = bridged_model(T; variant = var)
            x, y = MOI.add_variables(model, 2)
            c = MOI.add_constraint(model, vaf(T, [(1, x), (1, y)], -1), MOI.Nonnegatives(1))
            MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
            MOI.set(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}(), T(1) * x + T(1) * y)
            MOI.optimize!(model)
            @test MOI.get(model, MOI.TerminationStatus()) == MOI.OPTIMAL
            @test MOI.get(model, MOI.ObjectiveValue()) ≈ 1 atol = tol
            vx = MOI.get(model, MOI.VariablePrimal(), x)
            vy = MOI.get(model, MOI.VariablePrimal(), y)
            @test vx + vy ≈ 1 atol = tol
            @test iszero(vx) || iszero(vy) # one of them was fixed to zero
            @test MOI.get(model, MOI.ConstraintDual(), c) ≈ [1] atol = tol
        end

        @testset "dual infeasible" begin
            # minimize x - y subject to x + y ≥ 1 is unbounded; the presolve finds a ray
            model, optimizer = bridged_model(T; variant = var)
            x, y = MOI.add_variables(model, 2)
            MOI.add_constraint(model, vaf(T, [(1, x), (1, y)], -1), MOI.Nonnegatives(1))
            MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
            MOI.set(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}(), T(1) * x - T(1) * y)
            MOI.optimize!(model)
            @test optimizer.phasevalue == :presolve_dINF
            @test MOI.get(model, MOI.TerminationStatus()) == MOI.DUAL_INFEASIBLE
            @test MOI.get(model, MOI.PrimalStatus()) == MOI.INFEASIBILITY_CERTIFICATE
            @test MOI.get(model, MOI.DualStatus()) == MOI.NO_SOLUTION
            dx = MOI.get(model, MOI.VariablePrimal(), x)
            dy = MOI.get(model, MOI.VariablePrimal(), y)
            @test dx + dy ≈ 0 atol = 100 * eps(T)
            @test dx - dy < 0
            @test MOI.get(model, MOI.ObjectiveValue()) ≈ dx - dy
        end

        @testset "equality constraint with PSD constraint, sense $sense" for sense in (MOI.MIN_SENSE, MOI.MAX_SENSE)
            # minimize x1 + x3 subject to [x1 x2; x2 x3] ⪰ 0, x2 == 1.
            # The solution is x = (1, 1, 1). The dual of the PSD constraint is the
            # triangle (1, -1, 1) of the matrix [1 -1; -1 1], and the dual of the
            # equality constraint is 2: the objective is 2t for x2 == t.
            model, optimizer = bridged_model(T; variant = var)
            x = MOI.add_variables(model, 3)
            psd = MOI.add_constraint(model, MOI.VectorOfVariables(x), MOI.PositiveSemidefiniteConeTriangle(2))
            eq = MOI.add_constraint(model, T(1) * x[2], MOI.EqualTo(T(1)))
            s = sense == MOI.MIN_SENSE ? 1 : -1
            MOI.set(model, MOI.ObjectiveSense(), sense)
            MOI.set(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}(), T(s) * x[1] + T(s) * x[3] + T(5))
            MOI.optimize!(model)
            @test MOI.get(model, MOI.TerminationStatus()) == MOI.OPTIMAL
            @test MOI.get(model, MOI.PrimalStatus()) == MOI.FEASIBLE_POINT
            @test MOI.get(model, MOI.DualStatus()) == MOI.FEASIBLE_POINT
            @test MOI.get(model, MOI.ObjectiveValue()) ≈ 2s + 5 atol = tol
            @test MOI.get(model, MOI.DualObjectiveValue()) ≈ 2s + 5 atol = tol
            @test MOI.get(model, MOI.VariablePrimal(), x) ≈ [1, 1, 1] atol = tol
            @test MOI.get(model, MOI.ConstraintPrimal(), eq) ≈ 1 atol = tol
            @test MOI.get(model, MOI.ConstraintDual(), psd) ≈ [1, -1, 1] atol = tol
            @test MOI.get(model, MOI.ConstraintDual(), eq) ≈ 2 atol = tol
        end

        @testset "inconsistent equality constraints" begin
            model, optimizer = bridged_model(T; variant = var)
            x = MOI.add_variables(model, 2)
            MOI.add_constraint(model, MOI.VectorOfVariables(x), MOI.Nonnegatives(2))
            c1 = MOI.add_constraint(model, T(1) * x[1] + T(1) * x[2], MOI.EqualTo(T(1)))
            c2 = MOI.add_constraint(model, T(2) * x[1] + T(2) * x[2], MOI.EqualTo(T(3)))
            MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
            MOI.set(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}(), T(1) * x[1])
            MOI.optimize!(model)
            @test optimizer.phasevalue == :presolve_pINF
            @test MOI.get(model, MOI.TerminationStatus()) == MOI.INFEASIBLE
            @test MOI.get(model, MOI.PrimalStatus()) == MOI.NO_SOLUTION
            @test MOI.get(model, MOI.DualStatus()) == MOI.INFEASIBILITY_CERTIFICATE
            y1 = MOI.get(model, MOI.ConstraintDual(), c1)
            y2 = MOI.get(model, MOI.ConstraintDual(), c2)
            # `A' y = 0` and `b' y < 0`, where `b = (-1, -3)`
            @test y1 + 2y2 ≈ 0 atol = 100 * eps(T)
            @test -y1 - 3y2 < 0
            @test MOI.get(model, MOI.DualObjectiveValue()) ≈ y1 + 3y2
        end

        @testset "all variables fixed by equality constraints" begin
            for (rhs, status) in ((1, MOI.OPTIMAL), (-1, MOI.INFEASIBLE))
                model, optimizer = bridged_model(T; variant = var)
                x = MOI.add_variables(model, 2)
                MOI.add_constraint(model, MOI.VectorOfVariables(x), MOI.Nonnegatives(2))
                c1 = MOI.add_constraint(model, T(1) * x[1], MOI.EqualTo(T(rhs)))
                c2 = MOI.add_constraint(model, T(1) * x[1] + T(1) * x[2], MOI.EqualTo(T(3)))
                MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
                MOI.set(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}(), T(1) * x[1] + T(2) * x[2])
                MOI.optimize!(model)
                @test MOI.get(model, MOI.TerminationStatus()) == status
                @test MOI.get(model, MOI.VariablePrimal(), x) ≈ [rhs, 3 - rhs]
                if status == MOI.OPTIMAL
                    @test optimizer.phasevalue == :presolve_OPT
                    @test MOI.get(model, MOI.ObjectiveValue()) ≈ rhs + 2 * (3 - rhs)
                    @test MOI.get(model, MOI.DualStatus()) == MOI.FEASIBLE_POINT
                    # c = A_E' y with A_E = [1 0; 1 1] and c = (1, 2)
                    @test MOI.get(model, MOI.ConstraintDual(), c1) ≈ -1
                    @test MOI.get(model, MOI.ConstraintDual(), c2) ≈ 2
                else
                    # `x[2] = 4` and `x[1] = -1`: the constraint `x[1] >= 0` is violated,
                    # and its unit vector, mapped back to the equality constraints, is a
                    # certificate of infeasibility.
                    @test optimizer.phasevalue == :presolve_INF_cert
                    @test MOI.get(model, MOI.DualStatus()) == MOI.INFEASIBILITY_CERTIFICATE
                    @test MOI.get(model, MOI.DualObjectiveValue()) > 0
                end
            end
        end
    end

    # Regression tests for decisions that depend on rounding errors and on the scale of
    # the data. `T(1.3)` and `T(0.9)` are the `Float64` numbers converted to `T`, so
    # `1.3 * (0.9 / 1.3) - 0.9` is a rounding error, not zero.
    @testset "rounding errors and scaling, variant $var and type $T" for var in variants, T in (Float64, BigFloat)
        tol = soltol(var, T)
        f(T, terms, c = 0) = sum(T(a) * x for (a, x) in terms) + T(c)
        @testset "cancellation after substitution, zero cost" begin
            # min z s.t. 1.3x - 0.9y == 0, 1.3x - 0.9y <= 1, z >= 1: the solution has z = 1
            model, optimizer = bridged_model(T; variant = var)
            x, y, z = MOI.add_variables(model, 3)
            MOI.add_constraint(model, f(T, [(1.3, x), (-0.9, y)]), MOI.EqualTo(T(0)))
            MOI.add_constraint(model, f(T, [(1.3, x), (-0.9, y)]), MOI.LessThan(T(1)))
            MOI.add_constraint(model, T(1) * z, MOI.GreaterThan(T(1)))
            MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
            MOI.set(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}(), T(1) * z)
            MOI.optimize!(model)
            @test MOI.get(model, MOI.TerminationStatus()) == MOI.OPTIMAL
            @test MOI.get(model, MOI.ObjectiveValue()) ≈ 1 atol = tol
            vx, vy, vz = MOI.get(model, MOI.VariablePrimal(), [x, y, z])
            @test abs(T(1.3) * vx - T(0.9) * vy) <= tol
            @test abs(vx) + abs(vy) <= 1 # the variable fixed to zero is not left free
        end
        @testset "cancellation after substitution, unbounded" begin
            # min y s.t. 1.3x - 0.9y == 0, 1.3x - 0.9y <= 1 is unbounded
            model, optimizer = bridged_model(T; variant = var)
            x, y = MOI.add_variables(model, 2)
            MOI.add_constraint(model, f(T, [(1.3, x), (-0.9, y)]), MOI.EqualTo(T(0)))
            MOI.add_constraint(model, f(T, [(1.3, x), (-0.9, y)]), MOI.LessThan(T(1)))
            MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
            MOI.set(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}(), T(1) * y)
            MOI.optimize!(model)
            @test MOI.get(model, MOI.TerminationStatus()) == MOI.DUAL_INFEASIBLE
            @test MOI.get(model, MOI.PrimalStatus()) == MOI.INFEASIBILITY_CERTIFICATE
            dx, dy = MOI.get(model, MOI.VariablePrimal(), [x, y])
            # a ray: A_E d = 0, A_K d in the recession cone, c'd < 0
            @test abs(T(1.3) * dx - T(0.9) * dy) <= 100 * eps(T) * (abs(dx) + abs(dy))
            @test dy < 0
            @test MOI.get(model, MOI.ObjectiveValue()) ≈ dy
        end
        @testset "consistency of the objective with a large unrelated cost" begin
            # min 1e-6 x + 100 w s.t. x + y >= 1, w >= 1: x and y have equal columns and
            # different costs, so the problem is unbounded
            model, optimizer = bridged_model(T; variant = var)
            x, y, w = MOI.add_variables(model, 3)
            MOI.add_constraint(model, f(T, [(1, x), (1, y)]), MOI.GreaterThan(T(1)))
            MOI.add_constraint(model, T(1) * w, MOI.GreaterThan(T(1)))
            MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
            MOI.set(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}(), f(T, [(1e-6, x), (100, w)]))
            MOI.optimize!(model)
            @test optimizer.phasevalue == :presolve_dINF
            @test MOI.get(model, MOI.TerminationStatus()) == MOI.DUAL_INFEASIBLE
            dx, dy, dw = MOI.get(model, MOI.VariablePrimal(), [x, y, w])
            @test dx + dy == 0
            @test dw == 0
            @test T(1e-6) * dx < 0
        end
        @testset "consistency of equality constraints with a large unrelated right-hand side" begin
            # z == 1e6, x + y == 0, x + y == 1e-3 is inconsistent
            model, optimizer = bridged_model(T; variant = var)
            x, y, z = MOI.add_variables(model, 3)
            c1 = MOI.add_constraint(model, T(1) * z, MOI.EqualTo(T(1e6)))
            c2 = MOI.add_constraint(model, f(T, [(1, x), (1, y)]), MOI.EqualTo(T(0)))
            c3 = MOI.add_constraint(model, f(T, [(1, x), (1, y)]), MOI.EqualTo(T(1e-3)))
            MOI.add_constraint(model, MOI.VectorOfVariables([x, y]), MOI.Nonnegatives(2))
            MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
            MOI.set(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}(), f(T, [(1, x), (1, y)]))
            MOI.optimize!(model)
            @test optimizer.phasevalue == :presolve_pINF
            @test MOI.get(model, MOI.TerminationStatus()) == MOI.INFEASIBLE
            y1, y2, y3 = MOI.get(model, MOI.ConstraintDual(), [c1, c2, c3])
            # A'y = 0 and b'y < 0 with b = (-1e6, 0, -1e-3)
            @test y1 == 0
            @test abs(y2 + y3) <= 100 * eps(T) * (abs(y2) + abs(y3))
            @test -T(1e6) * y1 - T(1e-3) * y3 < 0
        end
        @testset "violated constraint when all variables are fixed" begin
            # x == 1 and [1000 - x; x - 1.000001] >= 0 is infeasible by 1e-6
            model, optimizer = bridged_model(T; variant = var)
            x = MOI.add_variable(model)
            ceq = MOI.add_constraint(model, T(1) * x, MOI.EqualTo(T(1)))
            cnn = MOI.add_constraint(model, MOI.Utilities.vectorize([T(1000) - T(1) * x, T(1) * x - parse(T, "1.000001")]), MOI.Nonnegatives(2))
            MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
            MOI.set(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}(), T(1) * x)
            MOI.optimize!(model)
            @test MOI.get(model, MOI.TerminationStatus()) == MOI.INFEASIBLE
            @test MOI.get(model, MOI.DualStatus()) == MOI.INFEASIBILITY_CERTIFICATE
            ye = MOI.get(model, MOI.ConstraintDual(), ceq)
            yn = MOI.get(model, MOI.ConstraintDual(), cnn)
            # A'y = 0, y in the dual cone, b'y < 0 with b = (-1, 1000, -1.000001)
            @test abs(ye - yn[1] + yn[2]) <= 100 * eps(T)
            @test all(yn .>= 0)
            @test -ye + 1000 * yn[1] - parse(T, "1.000001") * yn[2] < 0
        end
        @testset "primal ray from SDPA after the substitution" begin
            # min x + 3 s.t. x + y == 1, y >= 0 is unbounded; SDPA finds the ray for `y`
            model, optimizer = bridged_model(T; variant = var)
            x, y = MOI.add_variables(model, 2)
            MOI.add_constraint(model, f(T, [(1, x), (1, y)]), MOI.EqualTo(T(1)))
            MOI.add_constraint(model, T(1) * y, MOI.GreaterThan(T(0)))
            MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
            MOI.set(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}(), f(T, [(1, x)], 3))
            MOI.optimize!(model)
            @test MOI.get(model, MOI.TerminationStatus()) == MOI.DUAL_INFEASIBLE
            @test MOI.get(model, MOI.PrimalStatus()) == MOI.INFEASIBILITY_CERTIFICATE
            dx, dy = MOI.get(model, MOI.VariablePrimal(), [x, y])
            @test abs(dx + dy) <= 100 * eps(T) * (abs(dx) + abs(dy))
            @test dy >= 0
            @test dx < 0
            @test MOI.get(model, MOI.ObjectiveValue()) ≈ dx rtol = tol
        end
    end

    @testset "Convex: quantum state discrimination with variant $var" for var in (:sdpa_gmp,)
        ρ₁ = Complex{BigFloat}[1 0; 0 0]
        ρ₂ = (1//2) * Complex{BigFloat}[1 -im; im 1]
        E₁ = ComplexVariable(2, 2)
        E₂ = ComplexVariable(2, 2)
        problem = maximize(real((1//2) * tr(ρ₁ * E₁) + (1//2) * tr(ρ₂ * E₂)),
                           [E₁ ⪰ 0, E₂ ⪰ 0, E₁ + E₂ == Diagonal(ones(2))];
                           numeric_type = BigFloat)
        solve!(problem, () -> SDPAFamily.Optimizer(presolve = true, variant = var, silent = true))
        @test problem.optval ≈ 1//2 + 1 / (2 * sqrt(big(2))) atol = 1e-25
    end
end
