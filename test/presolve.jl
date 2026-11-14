using Test
using Convex
import MathOptInterface as MOI
using SparseArrays
using LinearAlgebra

import SDPAFamily
using Random
using Logging

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

include("presolve_kernel.jl")

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
        # Column 6 is a combination of columns 1, 2 and 5 up to rounding to `T`, which
        # counts as dependent (relative error at most `64 eps(T)`).
        A = sparse(T.(rand(-5:5, 40, 6)))
        A[:, 6] = A[:, 1] / 3 + A[:, 2] * T(7) / 10 - A[:, 5]
        c = transpose(A) * T.(rand(40))
        result = SDPAFamily.find_redundant_columns(A, c)
        @test result.redundant == [6]
        @test result.ray === nothing
        # a perturbation far above the rounding errors makes the columns independent
        A[1, 6] += T(1) / 10^6
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
            @test optimizer.presolve_data.mode in (:reduced, :none) # no fallback
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
            @test optimizer.presolve_data.mode in (:reduced, :none) # no fallback
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
            @test optimizer.presolve_data.mode in (:reduced, :none) # no fallback
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
            @test optimizer.presolve_data.mode in (:reduced, :none) # no fallback
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
                @test optimizer.presolve_data.mode in (:reduced, :none) # no fallback
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

    # Decisions that depend on rounding errors and on the scale of the data. `T(1.3)` and `T(0.9)` are the `Float64` numbers converted to `T`, so
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
            @test optimizer.presolve_data.mode in (:reduced, :none) # no fallback
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
            @test optimizer.presolve_data.mode in (:reduced, :none) # no fallback
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
            @test optimizer.presolve_data.mode in (:reduced, :none) # no fallback
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
            @test optimizer.presolve_data.mode in (:reduced, :none) # no fallback
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
            @test optimizer.presolve_data.mode in (:reduced, :none) # no fallback
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
            @test optimizer.presolve_data.mode in (:reduced, :none) # no fallback
            @test MOI.get(model, MOI.TerminationStatus()) == MOI.DUAL_INFEASIBLE
            @test MOI.get(model, MOI.PrimalStatus()) == MOI.INFEASIBILITY_CERTIFICATE
            dx, dy = MOI.get(model, MOI.VariablePrimal(), [x, y])
            @test abs(dx + dy) <= 100 * eps(T) * (abs(dx) + abs(dy))
            @test dy >= 0
            @test dx < 0
            @test MOI.get(model, MOI.ObjectiveValue()) ≈ dx rtol = tol
        end
    end

    @testset "exact dependence, zero blocks and small exact coefficients, variant $var and type $T" for var in variants, T in (Float64, BigFloat)
        tol = soltol(var, T)
        f(T, terms, c = 0) = sum(T(a) * x for (a, x) in terms) + T(c)
        @testset "exactly dependent equality row with a large multiplier" begin
            # 3x1 + 7x2 == 1, 100x1 + 11x2 + 13x3 == 2, row 3 = L * row 1 + row 2,
            # x >= 0, x3 <= 10, x3 <= 20, max x3: the solution is x = (0, 1/7, 3/91)
            L = T == Float64 ? T(10)^5 : T(10)^20
            model, optimizer = bridged_model(T; variant = var)
            x = MOI.add_variables(model, 3)
            MOI.add_constraint(model, f(T, [(3, x[1]), (7, x[2])]), MOI.EqualTo(T(1)))
            MOI.add_constraint(model, f(T, [(100, x[1]), (11, x[2]), (13, x[3])]), MOI.EqualTo(T(2)))
            MOI.add_constraint(model, f(T, [(3L + 100, x[1]), (7L + 11, x[2]), (13, x[3])]), MOI.EqualTo(L + 2))
            MOI.add_constraint(model, MOI.VectorOfVariables(x), MOI.Nonnegatives(3))
            MOI.add_constraint(model, T(1) * x[3], MOI.LessThan(T(10)))
            MOI.add_constraint(model, T(1) * x[3], MOI.LessThan(T(20)))
            MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
            MOI.set(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}(), T(1) * x[3])
            MOI.optimize!(model)
            @test optimizer.presolve_data.mode in (:reduced, :none) # no fallback
            @test MOI.get(model, MOI.TerminationStatus()) == MOI.OPTIMAL
            @test MOI.get(model, MOI.ObjectiveValue()) ≈ T(3) / 91 atol = tol
            @test MOI.get(model, MOI.VariablePrimal(), x) ≈ [0, T(1) / 7, T(3) / 91] atol = tol
        end
        @testset "zero PSD block at the determined point" begin
            # (a) X psd, X == 0; (b) X == [1 0; 0 0]; (c) constant zero PSD constraint;
            # (d) X == Y, Y == 0. All are feasible.
            for (name, rhs, chain) in (("a", [0, 0, 0], false), ("b", [1, 0, 0], false),
                                       ("d", [0, 0, 0], true))
                model, optimizer = bridged_model(T; variant = var)
                X = MOI.add_variables(model, 3)
                MOI.add_constraint(model, MOI.VectorOfVariables(X), MOI.PositiveSemidefiniteConeTriangle(2))
                if chain
                    Y = MOI.add_variables(model, 3)
                    for i in 1:3
                        MOI.add_constraint(model, T(1) * X[i] - T(1) * Y[i], MOI.EqualTo(T(0)))
                        MOI.add_constraint(model, T(1) * Y[i], MOI.EqualTo(T(0)))
                    end
                else
                    for i in 1:3
                        MOI.add_constraint(model, T(1) * X[i], MOI.EqualTo(T(rhs[i])))
                    end
                end
                MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
                MOI.set(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}(), T(1) * X[1])
                MOI.optimize!(model)
                @test optimizer.presolve_data.mode in (:reduced, :none) # no fallback
                @test MOI.get(model, MOI.TerminationStatus()) == MOI.OPTIMAL
                @test optimizer.phasevalue == :presolve_OPT
                @test MOI.get(model, MOI.ObjectiveValue()) == rhs[1]
            end
            model, optimizer = bridged_model(T; variant = var)
            x = MOI.add_variable(model)
            MOI.add_constraint(model, MOI.VectorAffineFunction(MOI.VectorAffineTerm{T}[], zeros(T, 3)),
                               MOI.PositiveSemidefiniteConeTriangle(2))
            MOI.add_constraint(model, T(1) * x, MOI.EqualTo(T(2)))
            MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
            MOI.set(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}(), T(1) * x)
            MOI.optimize!(model)
            @test optimizer.presolve_data.mode in (:reduced, :none) # no fallback
            @test MOI.get(model, MOI.TerminationStatus()) == MOI.OPTIMAL
            @test MOI.get(model, MOI.ObjectiveValue()) == 2
        end
        @testset "small coefficient after an exact cancellation" begin
            # x == y, w == y, L x - L y + s w + 1 >= 0, min y: the optimum is -1/s. The
            # coefficient s of the free variable must not be dropped.
            L, s = T == Float64 ? (T(1e4), T(1e-4)) : (big"1e40", big"1e-40")
            model, optimizer = bridged_model(T; variant = var)
            x, y, w = MOI.add_variables(model, 3)
            MOI.add_constraint(model, T(1) * x - T(1) * y, MOI.EqualTo(T(0)))
            MOI.add_constraint(model, T(1) * w - T(1) * y, MOI.EqualTo(T(0)))
            MOI.add_constraint(model, L * x - L * y + s * w, MOI.GreaterThan(T(-1)))
            MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
            MOI.set(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}(), T(1) * y)
            MOI.optimize!(model)
            @test optimizer.presolve_data.mode in (:reduced, :none) # no fallback
            # The presolve keeps `s`, and the scaling of the reduced problem
            # `min y : s y + 1 >= 0` makes it well scaled for all variants except
            # `:sdpa` with BigFloat data of size 1e±40.
            if !(var == :sdpa && T == BigFloat)
                @test MOI.get(model, MOI.TerminationStatus()) == MOI.OPTIMAL
                @test MOI.get(model, MOI.ObjectiveValue()) ≈ -1 / s rtol = sqrt(tol)
            end
        end
    end

    @testset "rows and columns dependent up to rounding of the data, $T" for T in (Float64, BigFloat)
        # 0.1x + 0.3y == 0.4 and 0.3x + 0.9y == 1.3 are dependent up to rounding to
        # `Float64` and inconsistent. With `Float64`, the presolve decides this with the
        # default tolerance; with `BigFloat` data that were rounded to `Float64`, with
        # `presolve_tol = 64 eps(Float64)`. The certificate is validated at the
        # tolerance of the decision, so the presolve returns it instead of falling back.
        ptol = T == Float64 ? nothing : 64 * eps(Float64)
        model, optimizer = bridged_model(T; variant = :sdpa_gmp, presolve_tol = ptol)
        x, y = MOI.add_variables(model, 2)
        c1 = MOI.add_constraint(model, T(0.1) * x + T(0.3) * y, MOI.EqualTo(T(0.4)))
        c2 = MOI.add_constraint(model, T(0.3) * x + T(0.9) * y, MOI.EqualTo(T(1.3)))
        MOI.add_constraint(model, MOI.VectorOfVariables([x, y]), MOI.Nonnegatives(2))
        MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
        MOI.set(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}(), T(1) * x)
        MOI.optimize!(model)
        @test optimizer.presolve_data.mode == :none
        @test MOI.get(model, MOI.TerminationStatus()) == MOI.INFEASIBLE
        @test MOI.get(model, MOI.DualStatus()) == MOI.INFEASIBILITY_CERTIFICATE
        # columns (0.1, 0.3) and (0.3, 0.9) with costs (1, 1): unbounded
        model, optimizer = bridged_model(T; variant = :sdpa_gmp, presolve_tol = ptol)
        x, y = MOI.add_variables(model, 2)
        MOI.add_constraint(model, T(0.1) * x + T(0.3) * y, MOI.GreaterThan(T(1)))
        MOI.add_constraint(model, T(0.3) * x + T(0.9) * y, MOI.GreaterThan(T(1)))
        MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
        MOI.set(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}(), T(1) * x + T(1) * y)
        MOI.optimize!(model)
        @test optimizer.presolve_data.mode == :none
        @test MOI.get(model, MOI.TerminationStatus()) == MOI.DUAL_INFEASIBLE
    end

    @testset "exact data that are Float64 numbers, BigFloat" begin
        # The data are exact; that they are also `Float64` numbers does not make them
        # rounded, so the default tolerance is 64 eps(BigFloat).
        T = BigFloat
        function solve(build!)
            model, optimizer = bridged_model(T; variant = :sdpa_gmp)
            build!(model)
            MOI.optimize!(model)
            return model, optimizer
        end
        # (a) (1e8+1)x + 1e8y == 1, 1e8x + (1e8-1)y == 0 (determinant -1), min x + y
        model, optimizer = solve(m -> begin
            x, y = MOI.add_variables(m, 2)
            MOI.add_constraint(m, T(10^8 + 1) * x + T(10^8) * y, MOI.EqualTo(T(1)))
            MOI.add_constraint(m, T(10^8) * x + T(10^8 - 1) * y, MOI.EqualTo(T(0)))
            MOI.add_constraint(m, T(1) * x + T(1) * y, MOI.GreaterThan(T(-10)))
            MOI.set(m, MOI.ObjectiveSense(), MOI.MIN_SENSE)
            MOI.set(m, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}(), T(1) * x + T(1) * y)
        end)
        @test optimizer.presolve_data.mode in (:reduced, :none) # no fallback
        @test MOI.get(model, MOI.TerminationStatus()) == MOI.OPTIMAL
        @test MOI.get(model, MOI.ObjectiveValue()) ≈ 1 atol = 1e-25
        # (b) x + y == 1 and x + y == 1 + 2^-50, x, y >= 0: infeasible
        model, optimizer = solve(m -> begin
            x, y = MOI.add_variables(m, 2)
            MOI.add_constraint(m, T(1) * x + T(1) * y, MOI.EqualTo(T(1)))
            MOI.add_constraint(m, T(1) * x + T(1) * y, MOI.EqualTo(1 + T(2)^-50))
            MOI.add_constraint(m, MOI.VectorOfVariables([x, y]), MOI.Nonnegatives(2))
            MOI.set(m, MOI.ObjectiveSense(), MOI.MIN_SENSE)
            MOI.set(m, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}(), T(1) * x)
        end)
        @test optimizer.presolve_data.mode in (:reduced, :none) # no fallback
        @test MOI.get(model, MOI.TerminationStatus()) == MOI.INFEASIBLE
        # (c) columns (1, 1) and (1, 1 + 2^-50) with costs (1, 1 + 2^-49): unbounded
        model, optimizer = solve(m -> begin
            x, y = MOI.add_variables(m, 2)
            MOI.add_constraint(m, T(1) * x + T(1) * y, MOI.GreaterThan(T(1)))
            MOI.add_constraint(m, T(1) * x + (1 + T(2)^-50) * y, MOI.GreaterThan(T(1)))
            MOI.set(m, MOI.ObjectiveSense(), MOI.MIN_SENSE)
            MOI.set(m, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}(), T(1) * x + (1 + T(2)^-49) * y)
        end)
        @test optimizer.presolve_data.mode in (:reduced, :none) # no fallback
        @test MOI.get(model, MOI.TerminationStatus()) == MOI.DUAL_INFEASIBLE
    end

    @testset "badly scaled small problems, $T and $var" for (T, var) in ((Float64, :sdpa), (Float64, :sdpa_gmp), (BigFloat, :sdpa_gmp))
        tol = soltol(var, T)
        # (A) x + y == 1, [K; x; y] >= 0 with a row without coefficients, min x + 2y:
        # the optimum is 1 at (1, 0) for every K >= 0.
        for K in (1e4, 1e10)
            model, optimizer = bridged_model(T; variant = var)
            x, y = MOI.add_variables(model, 2)
            MOI.add_constraint(model, T(1) * x + T(1) * y, MOI.EqualTo(T(1)))
            MOI.add_constraint(model, MOI.Utilities.operate(vcat, T, T(K), T(1) * x, T(1) * y), MOI.Nonnegatives(3))
            MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
            MOI.set(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}(), T(1) * x + T(2) * y)
            MOI.optimize!(model)
            @test optimizer.presolve_data.mode in (:reduced, :none) # no fallback
            @test MOI.get(model, MOI.TerminationStatus()) == MOI.OPTIMAL
            @test MOI.get(model, MOI.ObjectiveValue()) ≈ 1 atol = tol
            @test MOI.get(model, MOI.VariablePrimal(), [x, y]) ≈ [1, 0] atol = tol
        end
        # (B) x1 == x2, x1 >= R, x2 >= 0, min (x1 + x2) / R: the optimum is 2 at (R, R).
        for R in (1e10, 1e30)
            model, optimizer = bridged_model(T; variant = var)
            x = MOI.add_variables(model, 2)
            MOI.add_constraint(model, T(1) * x[1] - T(1) * x[2], MOI.EqualTo(T(0)))
            MOI.add_constraint(model, T(1) * x[1], MOI.GreaterThan(T(R)))
            MOI.add_constraint(model, T(1) * x[2], MOI.GreaterThan(T(0)))
            MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
            MOI.set(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}(), x[1] / T(R) + x[2] / T(R))
            MOI.optimize!(model)
            @test optimizer.presolve_data.mode in (:reduced, :none) # no fallback
            @test MOI.get(model, MOI.TerminationStatus()) == MOI.OPTIMAL
            @test MOI.get(model, MOI.ObjectiveValue()) ≈ 2 rtol = tol
            @test MOI.get(model, MOI.VariablePrimal(), x) ≈ [R, R] rtol = tol
        end
        # (C) s x1 + x2 == 1, x >= 0, min x2 - s x1: the optimum is -1 at (1/s, 0). The
        # coefficient s must not be dropped as rounding residue next to the 1.
        function caseC(s)
            model, optimizer = bridged_model(T; variant = var)
            x = MOI.add_variables(model, 2)
            MOI.add_constraint(model, T(s) * x[1] + T(1) * x[2], MOI.EqualTo(T(1)))
            MOI.add_constraint(model, MOI.VectorOfVariables(x), MOI.Nonnegatives(2))
            MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
            MOI.set(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}(), T(1) * x[2] - T(s) * x[1])
            MOI.optimize!(model)
            return model, optimizer
        end
        for s in (1e-20, 1e-300)
            model, optimizer = caseC(s)
            @test optimizer.presolve_data.mode in (:reduced, :none) # no fallback
            @test MOI.get(model, MOI.TerminationStatus()) == MOI.OPTIMAL
            @test MOI.get(model, MOI.ObjectiveValue()) ≈ -1 atol = tol
        end
        # With the subnormal s = 1e-310, the minimizer 1e310 is not a `Float64` number:
        # the result cannot be returned in `Float64`.
        if T == Float64
            model, optimizer = caseC(1e-310)
            @test optimizer.presolve_data.mode == :reduced # no fallback
            @test MOI.get(model, MOI.TerminationStatus()) == MOI.NUMERICAL_ERROR
            @test occursin("not finite", MOI.get(model, MOI.RawStatusString()))
        end
    end

    @testset "constant row with a tiny negative constant" begin
        # From Convex's `constant_fix!_with_vectors`: `y` is fixed to the value of a
        # previous solve, 0.5 - 2.9e-32, so `y >= 0.5` becomes a row without
        # coefficients and a tiny negative constant. Scaling that row by its constant
        # made SDPA-GMP report `pdINF`; rows without coefficients are not scaled.
        T = BigFloat
        x = ComplexVariable(5)
        fix!(x, ones(5) + im * ones(5))
        y = Variable()
        prob = minimize(real(y * sum(x)), [y >= 0.5, real(x) >= 0.5, imag(x) >= 0]; numeric_type = T)
        solve!(prob, () -> SDPAFamily.Optimizer{T}(silent = true, variant = :sdpa_gmp))
        free!(x)
        fix!(y)
        optimizer = SDPAFamily.Optimizer{T}(silent = true, variant = :sdpa_gmp)
        solve!(prob, () -> optimizer)
        @test optimizer.presolve_data.mode in (:reduced, :none) # no fallback
        @test optimizer.phasevalue == :pdOPT
        @test prob.optval ≈ 1.25 atol = 1e-20
    end

    @testset "nothing to reduce, variant $var, sense $sense" for var in (:sdpa_gmp, :sdpa_dd), sense in (MOI.MIN_SENSE, MOI.MAX_SENSE)
        # no equality constraints and independent columns: the model is passed unchanged
        T = BigFloat
        s = sense == MOI.MIN_SENSE ? 1 : -1
        model, optimizer = bridged_model(T; variant = var)
        x, y = MOI.add_variables(model, 2)
        MOI.add_constraint(model, MOI.VectorOfVariables([x, y]), MOI.Nonnegatives(2))
        MOI.add_constraint(model, T(1) * x + T(2) * y, MOI.GreaterThan(T(2)))
        MOI.set(model, MOI.ObjectiveSense(), sense)
        MOI.set(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}(), T(s) * x + T(s) * y)
        MOI.optimize!(model)
        @test optimizer.presolve_data.mode in (:reduced, :none) # no fallback
        data = optimizer.presolve_data
        @test data.mode == :reduced
        @test data.elim === nothing
        @test isempty(data.redundant)
        @test MOI.get(model, MOI.TerminationStatus()) == MOI.OPTIMAL
        @test MOI.get(model, MOI.ObjectiveValue()) ≈ s atol = 1e-10
    end

    @testset "fallback after the solve keeps the objective constant, variant $var" for var in (:sdpa_gmp, :sdpa_dd)
        # x + y == 3 with objective 2x + y: the substitution adds a constant to the
        # objective; a fallback after the solve must remove it again.
        T = BigFloat
        model, optimizer = bridged_model(T; variant = var)
        x, y = MOI.add_variables(model, 2)
        MOI.add_constraint(model, T(1) * x + T(1) * y, MOI.EqualTo(T(3)))
        MOI.add_constraint(model, MOI.VectorOfVariables([x, y]), MOI.Nonnegatives(2))
        MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
        MOI.set(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}(), T(2) * x + T(1) * y)
        MOI.optimize!(model)
        @test optimizer.presolve_data.mode in (:reduced, :none) # no fallback
        @test MOI.get(model, MOI.ObjectiveValue()) ≈ 3 atol = 1e-10
        # the same model through the fallback path, starting from the reduced model
        data = optimizer.presolve_data
        @test data.mode == :reduced
        # an intentional fallback, not counted by the logger of the tests
        fb = Logging.with_logger(Logging.NullLogger()) do
            SDPAFamily.fallback_model!(optimizer, "test")
        end
        SDPAFamily.run_sdpa!(optimizer, fb, joinpath(optimizer.tempdir, "input.dat-s"),
                             joinpath(optimizer.tempdir, "output.dat"))
        SDPAFamily.postsolve!(optimizer)
        @test data.mode == :fallback
        @test optimizer.objconstant == 0
        if MOI.get(optimizer, MOI.TerminationStatus()) == MOI.OPTIMAL
            @test MOI.get(optimizer, MOI.ObjectiveValue()) ≈ 3 atol = 1e-6
        end
    end

    @testset "dense equality constraints, $T" for (T, var, m) in ((Float64, :sdpa, 60), (BigFloat, :sdpa_gmp, 60))
        # m dense random equality constraints on m + 10 variables, x >= x0 - 1, min sum(x)
        rng = Random.MersenneTwister(5)
        n = m + 10
        A = T.(randn(rng, m, n))
        x0 = T.(randn(rng, n))
        model, optimizer = bridged_model(T; variant = var)
        x = MOI.add_variables(model, n)
        for i in 1:m
            MOI.add_constraint(model, MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(A[i, :], x), zero(T)),
                               MOI.EqualTo(sum(A[i, :] .* x0)))
        end
        for j in 1:n
            MOI.add_constraint(model, T(1) * x[j], MOI.GreaterThan(x0[j] - 1))
        end
        MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
        MOI.set(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}(),
                MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(one(T), x), zero(T)))
        MOI.optimize!(model)
        @test optimizer.presolve_data.mode in (:reduced, :none) # no fallback
        @test MOI.get(model, MOI.TerminationStatus()) == MOI.OPTIMAL
        v = MOI.get(model, MOI.VariablePrimal(), x)
        residual = abs.(A * v - A * x0) ./ (abs.(A) * abs.(v) + abs.(A) * abs.(x0))
        @test maximum(residual) <= 1000 * eps(T)
        @test all(v .>= x0 .- 1 .- soltol(var, T))
    end

    @testset "objective consistent with a free direction up to rounding, $T" for T in (Float64, BigFloat)
        # x == y, 0 <= z <= 1, min a x - b y + z with a = 0.1 + 0.2 and b = 0.3 (as
        # `Float64` numbers, a - b = 2^-54). The direction (1, 1, 0) changes no
        # constraint. With `Float64`, a - b is zero up to the data tolerance relative to
        # a and b, so the free variable is dropped. With `BigFloat`, the data are exact and
        # the problem is unbounded.
        model, optimizer = bridged_model(T; variant = :sdpa_gmp)
        x, y, z = MOI.add_variables(model, 3)
        MOI.add_constraint(model, T(1) * x - T(1) * y, MOI.EqualTo(T(0)))
        MOI.add_constraint(model, T(1) * z, MOI.GreaterThan(T(0)))
        MOI.add_constraint(model, T(1) * z, MOI.LessThan(T(1)))
        MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
        MOI.set(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}(),
                T(0.1 + 0.2) * x - T(0.3) * y + T(1) * z)
        MOI.optimize!(model)
        @test optimizer.presolve_data.mode in (:reduced, :none) # no fallback
        if T == Float64
            @test MOI.get(model, MOI.TerminationStatus()) == MOI.OPTIMAL
            @test MOI.get(model, MOI.ObjectiveValue()) ≈ 0 atol = 1e-6
        else
            @test MOI.get(model, MOI.TerminationStatus()) == MOI.DUAL_INFEASIBLE
            d = MOI.get(model, MOI.VariablePrimal(), [x, y, z])
            @test d[1] == d[2] && d[3] == 0 && d[1] < 0
        end
    end

    @testset "two parts at different scales, $T and $var" for (T, var) in ((Float64, :sdpa), (BigFloat, :sdpa_gmp))
        # min -2^-18 x + 2^-39 y with x <= 5 * 2^18 and 2^-10 y >= -3 * 2^30: the optimum
        # is -5 - 3 * 2^1 = -11 at (5 * 2^18, -3 * 2^40). With one factor for the whole
        # right-hand side, the first part has size 2^-12 in the scaled problem, below
        # the tolerance of `:sdpa`; each part gets its own factor.
        model, optimizer = bridged_model(T; variant = var)
        x, y = MOI.add_variables(model, 2)
        MOI.add_constraint(model, T(1) * x, MOI.LessThan(T(5 * 2^18)))
        MOI.add_constraint(model, T(2)^-10 * y, MOI.GreaterThan(-T(3) * T(2)^30))
        MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
        MOI.set(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}(), -T(2)^-18 * x + T(2)^-39 * y)
        MOI.optimize!(model)
        @test optimizer.presolve_data.mode in (:reduced, :none) # no fallback
        @test MOI.get(model, MOI.TerminationStatus()) == MOI.OPTIMAL
        @test MOI.get(model, MOI.ObjectiveValue()) ≈ -11 rtol = 10 * soltol(var, T)
        @test MOI.get(model, MOI.VariablePrimal(), [x, y]) ≈ [5 * 2^18, -3 * T(2)^40] rtol = 10 * soltol(var, T)
    end

    @testset "certificates from SDPA with many rows, $T and $var" for (T, var) in ((Float64, :sdpa), (BigFloat, :sdpa_gmp))
        # The presolve finds no reduction; SDPA reports the infeasibility. Its
        # certificates are iterates at its objective bound, and must be accepted
        # whatever the number of rows and variables.
        n = 120
        # x_i >= 1 for all i and sum(x) <= n / 2: infeasible
        model, optimizer = bridged_model(T; variant = var)
        x = MOI.add_variables(model, n)
        MOI.add_constraint(model, MOI.Utilities.vectorize([T(1) * xi - T(1) for xi in x]), MOI.Nonnegatives(n))
        MOI.add_constraint(model, sum(T(1) * xi for xi in x), MOI.LessThan(T(n) / 2))
        MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
        MOI.set(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}(), sum(T(1) * xi for xi in x))
        MOI.optimize!(model)
        @test optimizer.presolve_data.mode == :reduced # SDPA was called, no fallback
        @test MOI.get(model, MOI.TerminationStatus()) == MOI.INFEASIBLE
        # x_i >= 0 for all i and x_1 >= x_2, min -sum(x): unbounded
        model, optimizer = bridged_model(T; variant = var)
        x = MOI.add_variables(model, n)
        MOI.add_constraint(model, MOI.VectorOfVariables(x), MOI.Nonnegatives(n))
        MOI.add_constraint(model, T(1) * x[1] - T(1) * x[2], MOI.GreaterThan(T(0)))
        MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
        MOI.set(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}(), -sum(T(1) * xi for xi in x))
        MOI.optimize!(model)
        @test optimizer.presolve_data.mode == :reduced
        @test MOI.get(model, MOI.TerminationStatus()) == MOI.DUAL_INFEASIBLE
    end

    @testset "inactive variable with a large range, $T" for T in (Float64, BigFloat)
        # 1e150 x1 + x2 == 1e150, x >= 0, 1 <= x3 <= 5, min x1 + x2 + 1e-150 x3: the
        # optimum is 1 + 1e-150 at (1, 0, 1), but x2 ranges up to 1e150 on the feasible
        # set. With the first scaling the objective is 1e-150 in the scaled problem and
        # SDPA's gap tolerance does not resolve it; that result is rejected, and the
        # second scaling gives the solution.
        var = T == Float64 ? :sdpa : :sdpa_gmp
        model, optimizer = bridged_model(T; variant = var)
        x = MOI.add_variables(model, 3)
        MOI.add_constraint(model, T(1e150) * x[1] + T(1) * x[2], MOI.EqualTo(T(1e150)))
        MOI.add_constraint(model, MOI.VectorOfVariables(x[1:2]), MOI.Nonnegatives(2))
        MOI.add_constraint(model, T(1e-150) * x[3], MOI.GreaterThan(T(1e-150)))
        MOI.add_constraint(model, T(1) * x[3], MOI.LessThan(T(5)))
        MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
        MOI.set(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}(), T(1) * x[1] + T(1) * x[2] + T(1e-150) * x[3])
        MOI.optimize!(model)
        @test optimizer.presolve_data.mode == :reduced # no fallback
        @test optimizer.presolve_data.border # second scaling
        @test MOI.get(model, MOI.TerminationStatus()) == MOI.OPTIMAL
        @test MOI.get(model, MOI.ObjectiveValue()) ≈ 1 rtol = 10 * soltol(var, T)
    end

    @testset "bound that the first scaling makes tiny, $T and $var" for (T, var, ss) in
            ((Float64, :sdpa, (1e-4, 1e-12, 1e-16)), (BigFloat, :sdpa_gmp, (big"1e-10", big"1e-35", big"1e-50")))
        # s x1 + x2 == 1, 0 <= x1 <= 1, x2 >= 0, min x2 - s x1: the optimum is at x1 = 1.
        # Equilibrating the reduced matrix alone balances the coefficient s of
        # x2 = 1 - s x1 >= 0 against the bounds of x1, and the right-hand side of
        # x1 <= 1 falls below SDPA's tolerance in the scaled problem for small s. SDPA
        # then returns x1 far above 1. The check of the model rows rejects that point,
        # and the second scaling gives the solution.
        for s in ss
            model, optimizer = bridged_model(T; variant = var)
            x = MOI.add_variables(model, 2)
            MOI.add_constraint(model, T(s) * x[1] + T(1) * x[2], MOI.EqualTo(T(1)))
            MOI.add_constraint(model, T(1) * x[1], MOI.GreaterThan(T(0)))
            MOI.add_constraint(model, T(1) * x[1], MOI.LessThan(T(1)))
            MOI.add_constraint(model, T(1) * x[2], MOI.GreaterThan(T(0)))
            MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
            MOI.set(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}(), T(1) * x[2] - T(s) * x[1])
            MOI.optimize!(model)
            @test optimizer.presolve_data.mode == :reduced # no fallback
            @test MOI.get(model, MOI.TerminationStatus()) == MOI.OPTIMAL
            @test MOI.get(model, MOI.VariablePrimal(), x[1]) ≈ 1 atol = 1e-2
            @test MOI.get(model, MOI.ObjectiveValue()) ≈ 1 - 2s atol = 10 * soltol(var, T)
        end
    end

    @testset "bound far from the scale that both scalings give, variant $var" for var in (:sdpa, :sdpa_gmp)
        # The same model with 0 <= x1 <= 1e190 and s = 1e-200 (optimum at x1 = 1e190,
        # objective 1 - 2e-10). Both scalings give x1 a unit far below 1e190. The first
        # result violates the bound and is rejected; with the second, SDPA reports
        # that the problem is unbounded, which is not accepted after feasible iterates.
        T = Float64
        model, optimizer = bridged_model(T; variant = var)
        x = MOI.add_variables(model, 2)
        MOI.add_constraint(model, T(1e-200) * x[1] + T(1) * x[2], MOI.EqualTo(T(1)))
        MOI.add_constraint(model, T(1) * x[1], MOI.GreaterThan(T(0)))
        MOI.add_constraint(model, T(1) * x[1], MOI.LessThan(T(1e190)))
        MOI.add_constraint(model, T(1) * x[2], MOI.GreaterThan(T(0)))
        MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
        MOI.set(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}(), T(1) * x[2] - T(1e-200) * x[1])
        MOI.optimize!(model)
        @test MOI.get(model, MOI.TerminationStatus()) == MOI.NUMERICAL_ERROR
        @test MOI.get(model, MOI.PrimalStatus()) == MOI.UNKNOWN_RESULT_STATUS
    end

    @testset "optimum beyond SDPA's objective bounds" begin
        # x1 >= 1, x_{i+1} >= 2 x_i, z == x1, min sum(x): the optimum is 2^n - 1. SDPA
        # reports unboundedness when an objective value passes its parameters
        # `lowerBound`/`upperBound` (±1e5 by default), so with the defaults `:sdpa_dd`
        # reports this feasible problem as infeasible for n = 20. With larger bounds
        # it is solved.
        T, n = BigFloat, 20
        model, optimizer = bridged_model(T; variant = :sdpa_dd, params = (lowerBound = -1e10, upperBound = 1e10))
        x = MOI.add_variables(model, n)
        z = MOI.add_variable(model)
        MOI.add_constraint(model, T(1) * z - T(1) * x[1], MOI.EqualTo(T(0)))
        MOI.add_constraint(model, T(1) * x[1], MOI.GreaterThan(T(1)))
        for i in 1:n-1
            MOI.add_constraint(model, T(1) * x[i+1] - T(2) * x[i], MOI.GreaterThan(T(0)))
        end
        MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
        MOI.set(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}(), sum(T(1) * xi for xi in x))
        MOI.optimize!(model)
        @test MOI.get(model, MOI.TerminationStatus()) == MOI.OPTIMAL
        @test MOI.get(model, MOI.ObjectiveValue()) ≈ 2^n - 1 rtol = 1e-10
    end

    @testset "iteration limit with feasible iterates, sense $sense" for sense in (MOI.MIN_SENSE, MOI.MAX_SENSE)
        # SDPA stops after 4 iterations with the status `pdFEAS` and a relative gap of
        # about 1e-2. This is reported as `ITERATION_LIMIT`, not as `OPTIMAL`.
        T = Float64
        s = sense == MOI.MIN_SENSE ? 1 : -1
        model, optimizer = bridged_model(T; variant = :sdpa, params = (maxIteration = 4,))
        x, y, z = MOI.add_variables(model, 3)
        MOI.add_constraint(model, T(1) * x, MOI.GreaterThan(T(1)))
        MOI.add_constraint(model, T(1) * y, MOI.GreaterThan(T(1)))
        MOI.add_constraint(model, T(1) * x + T(1) * y + T(1) * z, MOI.EqualTo(T(3)))
        MOI.add_constraint(model, MOI.Utilities.operate(vcat, T, T(1) * z, T(0.5), T(1) * y),
                           MOI.PositiveSemidefiniteConeTriangle(2))
        MOI.set(model, MOI.ObjectiveSense(), sense)
        MOI.set(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}(), T(2s) * x + T(s) * y + T(3s) * z)
        MOI.optimize!(model)
        @test optimizer.presolve_data.mode == :reduced # no fallback
        @test optimizer.phasevalue == :pdFEAS_gap
        @test MOI.get(model, MOI.TerminationStatus()) == MOI.ITERATION_LIMIT
        @test MOI.get(model, MOI.PrimalStatus()) == MOI.FEASIBLE_POINT
        @test MOI.get(model, MOI.DualStatus()) == MOI.FEASIBLE_POINT
        v = MOI.get(model, MOI.VariablePrimal(), [x, y, z])
        @test v[1] + v[2] + v[3] ≈ 3 atol = 1e-6
        @test all(v[1:2] .>= 1 - 1e-6)
    end

    @testset "iteration limit with feasible iterates without presolve" begin
        # min 2x + y with x >= 1, y >= 1, stopped after 4 iterations: SDPA-GMP reports
        # `pdFEAS` with objective values 1.4e5 and 3. The gap criterion is the same
        # with and without the presolve.
        T = BigFloat
        model, optimizer = bridged_model(T; variant = :sdpa_gmp, presolve = false, params = (maxIteration = 4,))
        x, y = MOI.add_variables(model, 2)
        MOI.add_constraint(model, MOI.Utilities.operate(vcat, T, T(1) * x - T(1), T(1) * y - T(1)), MOI.Nonnegatives(2))
        MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
        MOI.set(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}(), T(2) * x + T(1) * y)
        MOI.optimize!(model)
        @test optimizer.phasevalue == :pdFEAS_gap
        @test MOI.get(model, MOI.TerminationStatus()) == MOI.ITERATION_LIMIT
        @test MOI.get(model, MOI.PrimalStatus()) == MOI.FEASIBLE_POINT
    end

    include("presolve_random.jl")
    @testset "random badly scaled problems, $T and $var" for (T, var) in ((BigFloat, :sdpa_gmp), (Float64, :sdpa))
        # Rows scaled by 10^±12 or 10^±6 and variables by 2^±20. The result must be
        # that of the same problem without the scaling of the variables.
        rng = Random.MersenneTwister(12)
        tol = 100 * soltol(var, T)
        for k in 1:20
            p = random_conic_problem(rng, T; scale = true)
            s = T(2) .^ rand(rng, (-20, 0, 20), length(p.c))
            objective = map((p, scale_variables(p, s))) do q
                model, optimizer = bridged_model(T; variant = var)
                add_random_conic_problem!(model, q, T)
                MOI.optimize!(model)
                @test optimizer.presolve_data.mode in (:reduced, :none) # no fallback
                @test MOI.get(model, MOI.TerminationStatus()) == MOI.OPTIMAL
                pobj, dobj = MOI.get(model, MOI.ObjectiveValue()), MOI.get(model, MOI.DualObjectiveValue())
                @test abs(pobj - dobj) <= tol * max(1, abs(pobj))
                pobj
            end
            @test objective[1] ≈ objective[2] atol = tol * max(1, abs(objective[1]))
        end
    end

    @testset "random infeasible and unbounded problems with scaled variables" begin
        # The status must be the one the problem was built for.
        T, var = BigFloat, :sdpa_gmp
        rng = Random.MersenneTwister(11)
        for k in 1:30
            kind = rand(rng, (:infeasible, :unbounded))
            p = random_conic_problem(rng, T; scale = true, kind = kind)
            p = scale_variables(p, T(2) .^ rand(rng, (-20, 0, 20), length(p.c)))
            model, optimizer = bridged_model(T; variant = var)
            add_random_conic_problem!(model, p, T)
            MOI.optimize!(model)
            @test optimizer.presolve_data.mode in (:reduced, :none) # no fallback
            @test MOI.get(model, MOI.TerminationStatus()) ==
                  (kind == :infeasible ? MOI.INFEASIBLE : MOI.DUAL_INFEASIBLE)
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
