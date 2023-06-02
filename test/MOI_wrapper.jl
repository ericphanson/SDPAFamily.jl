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
        # The output file is possibly corrupted.
        r"test_attribute_RawStatusString$",
        r"test_attribute_SolveTimeSec$",
        r"test_conic_empty_matrix$",
        r"test_solve_TerminationStatus_DUAL_INFEASIBLE$",
        # Unable to bridge RotatedSecondOrderCone to PSD because the dimension is too small: got 2, expected >= 3.
        r"test_conic_SecondOrderCone_INFEASIBLE$",
        r"test_constraint_PrimalStart_DualStart_SecondOrderCone$",
        # Expression: MOI.get(model, MOI.TerminationStatus()) == MOI.INFEASIBLE
        #  Evaluated: MathOptInterface.INFEASIBLE_OR_UNBOUNDED == MathOptInterface.INFEASIBLE
        r"test_conic_NormInfinityCone_INFEASIBLE$",
        r"test_conic_NormOneCone_INFEASIBLE$",
        r"test_solve_DualStatus_INFEASIBILITY_CERTIFICATE_Interval_lower$",
        r"test_solve_DualStatus_INFEASIBILITY_CERTIFICATE_Interval_upper$",
        # Incorrect objective
        # See https://github.com/jump-dev/MathOptInterface.jl/issues/1759
        r"test_unbounded_MIN_SENSE$",
        r"test_unbounded_MIN_SENSE_offset$",
        r"test_unbounded_MAX_SENSE$",
        r"test_unbounded_MAX_SENSE_offset$",
        r"test_infeasible_MAX_SENSE$",
        r"test_infeasible_MAX_SENSE_offset$",
        r"test_infeasible_MIN_SENSE$",
        r"test_infeasible_MIN_SENSE_offset$",
        r"test_infeasible_affine_MAX_SENSE$",
        r"test_infeasible_affine_MAX_SENSE_offset$",
        r"test_infeasible_affine_MIN_SENSE$",
        r"test_infeasible_affine_MIN_SENSE_offset$",
        # TODO remove when PR merged
        # See https://github.com/jump-dev/MathOptInterface.jl/pull/1769
        r"test_objective_ObjectiveFunction_blank$",
        # FIXME investigate
        #  Expression: isapprox(MOI.get(model, MOI.ObjectiveValue()), T(2), config)
        #   Evaluated: isapprox(5.999999984012059, 2.0, ...
        r"test_modification_delete_variables_in_a_batch$",
        # FIXME investigate
        #  Expression: isapprox(MOI.get(model, MOI.ObjectiveValue()), objective_value, config)
        #   Evaluated: isapprox(-2.1881334077988868e-7, 5.0, ...
        r"test_objective_qp_ObjectiveFunction_edge_case$",
        # FIXME investigate
        #  Expression: isapprox(MOI.get(model, MOI.ObjectiveValue()), objective_value, config)
        #   Evaluated: isapprox(-2.1881334077988868e-7, 5.0, ...
        r"test_objective_qp_ObjectiveFunction_zero_ofdiag$",
        # FIXME investigate
        #  Expression: isapprox(MOI.get(model, MOI.ConstraintPrimal(), index), solution_value, config)
        #   Evaluated: isapprox(2.5058846553349667e-8, 1.0, ...
        r"test_variable_solve_with_lowerbound$",
        # FIXME investigate
        # See https://github.com/jump-dev/SDPA.jl/runs/7246518765?check_suite_focus=true#step:6:128
        # Expression: ≈(MOI.get(model, MOI.ConstraintDual(), c), T[1, 0, 0, -1, 1, 0, -1, -1, 1] / T(3), config)
        #  Evaluated: ≈([0.3333333625488728, -0.16666659692134123, -0.16666659693012292, -0.16666659692134123, 0.33333336253987234, -0.16666659692112254, -0.16666659693012292, -0.16666659692112254, 0.333333362548654], [0.3333333333333333, 0.0, 0.0, -0.3333333333333333, 0.3333333333333333, 0.0, -0.3333333333333333, -0.3333333333333333, 0.3333333333333333]
        r"test_conic_PositiveSemidefiniteConeSquare_3$",
        # FIXME
        r"test_model_LowerBoundAlreadySet$",
        r"test_model_UpperBoundAlreadySet$",
    ]
    if var != :sdpa
        append!(
            exclude,
            [
                r"test_solve_VariableIndex_ConstraintDual_MAX_SENSE$",
                r"test_solve_VariableIndex_ConstraintDual_MIN_SENSE$",
                r"test_constraint_ScalarAffineFunction_EqualTo$",
                # ITERATION_LIMIT for sdpa_dd
                r"test_conic_SecondOrderCone_negative_post_bound_2$",
                r"test_conic_SecondOrderCone_negative_post_bound_3$",
                r"test_conic_SecondOrderCone_no_initial_bound$",
                r"test_linear_FEASIBILITY_SENSE$",
                r"test_linear_transform$",
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
            # BigFloat # not yet supported: MathOptInterface#41
        )
            MOI_tests(var, Float64)
            MOI_tests(var, BigFloat)
        end
    end
end

MOI_tests()
