using Test

import MathOptInterface as MOI

import SDPAFamily

function MOI_tests(var, ::Type{T}) where {T}
    optimizer = SDPAFamily.Optimizer{T}(presolve = true, variant = var)
    MOI.set(optimizer, MOI.Silent(), true)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("maxIteration"), 5000)
    @testset "SolverName" begin
        @test MOI.get(optimizer, MOI.SolverName()) == "SDPAFamily"
    end

    # UniversalFallback is needed for starting values, even if they are ignored by SDPA
    cache = MOI.Utilities.UniversalFallback(MOI.Utilities.Model{T}())
    cached = MOI.Utilities.CachingOptimizer(cache, optimizer)
    bridged = MOI.Bridges.full_bridge_optimizer(cached, T)
    # test 1e-3 because of rsoc3 test, otherwise, 1e-5 is enough
    config = MOI.Test.Config(
        T,
        atol=1e-3,
        rtol=1e-3,
        exclude = Any[
            MOI.ConstraintBasisStatus,
            MOI.VariableBasisStatus,
            MOI.ObjectiveBound,
            MOI.SolverVersion,
        ],
    )

    exclude = [
        # The infeasibility certificates returned from SDPA's output do not satisfy
        # the relations between ray, objective and duals that these tests check. See https://github.com/jump-dev/MathOptInterface.jl/issues/1759
        # Expression: isapprox(T(11 // 5) * d, obj, config)
        #  Evaluated: isapprox(-158862.07316834212, -10.999999999999995, ...
        r"test_unbounded_MIN_SENSE$",
        r"test_unbounded_MIN_SENSE_offset$",
        r"test_unbounded_MAX_SENSE$",
        r"test_unbounded_MAX_SENSE_offset$",
        # Expression: isapprox(dl + du, T(0), config)
        #  Evaluated: isapprox(-2.199999999254942, 0.0, ...
        r"test_infeasible_MAX_SENSE$",
        r"test_infeasible_MAX_SENSE_offset$",
        r"test_infeasible_MIN_SENSE$",
        r"test_infeasible_MIN_SENSE_offset$",
        r"test_infeasible_affine_MAX_SENSE$",
        r"test_infeasible_affine_MAX_SENSE_offset$",
        r"test_infeasible_affine_MIN_SENSE$",
        r"test_infeasible_affine_MIN_SENSE_offset$",
    ]
    if var == :sdpa
        append!(
            exclude,
            [
                # Expression: MOI.get(model, MOI.PrimalStatus()) in (MOI.NO_SOLUTION, MOI.INFEASIBLE_POINT)
                #  Evaluated: MathOptInterface.UNKNOWN_RESULT_STATUS in ...
                r"test_conic_SecondOrderCone_INFEASIBLE$",
                # Expression: MOI.get(model, MOI.TerminationStatus()) == config.infeasible_status
                #  Evaluated: MathOptInterface.INFEASIBLE_OR_UNBOUNDED == MathOptInterface.INFEASIBLE
                r"test_solve_DualStatus_INFEASIBILITY_CERTIFICATE_Interval_lower$",
                r"test_solve_DualStatus_INFEASIBILITY_CERTIFICATE_Interval_upper$",
            ],
        )
    end
    if var in (:sdpa_dd, :sdpa_qd)
        append!(
            exclude,
            [
                # The constraints `x >= 1` and `x <= 1` leave no strictly feasible point,
                # and the presolve only eliminates `EqualTo` and `Zeros` constraints.
                # Expression: ≈(sl + su, T(1), atol = config.atol)
                #  Evaluated: 0.0 ≈ 1.0 (atol=0.001)
                r"test_solve_VariableIndex_ConstraintDual_MAX_SENSE$",
                r"test_solve_VariableIndex_ConstraintDual_MIN_SENSE$",
            ],
        )
    end

    MOI.Test.runtests(
        bridged,
        config;
        exclude,
    )
end

function MOI_tests()
    @testset "MOI tests with variant $var" for var in variants
        @testset "MOI tests with type T=$T" for T in (
            Float64,
            BigFloat
        )
            MOI_tests(var, T)
        end
    end
end

MOI_tests()
