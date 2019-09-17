using Test

using MathOptInterface
const MOI = MathOptInterface
const MOIT = MOI.Test
const MOIU = MOI.Utilities
const MOIB = MOI.Bridges

import SDPAFamily

@testset "MOI tests" for var in variants
    @testset "MOI tests with variant $var and type T=$T" for T in (
                                                    Float64,
                                                    # BigFloat # not yet supported: MathOptInterface#41
                                                )
        optimizer = SDPAFamily.Optimizer{T}(presolve=true, variant = var)
        MOI.set(optimizer, MOI.Silent(), true)

        @testset "SolverName" begin
            @test MOI.get(optimizer, MOI.SolverName()) == "SDPAFamily"
        end

        @testset "supports_default_copy_to" begin
            @test MOIU.supports_allocate_load(optimizer, false)
            @test !MOIU.supports_allocate_load(optimizer, true)
        end

        # UniversalFallback is needed for starting values, even if they are ignored by SDPA
        cache = MOIU.UniversalFallback(MOIU.Model{T}())
        cached = MOIU.CachingOptimizer(cache, optimizer)
        bridged = MOIB.full_bridge_optimizer(cached, T)
        # test 1e-3 because of rsoc3 test, otherwise, 1e-5 is enough
        config = MOIT.TestConfig(atol=1e-3, rtol=1e-3)

        @testset "Unit" begin
            exclusion_list = [
                # `TimeLimitSec` not supported.
                "time_limit_sec",
                # SingleVariable objective of bridged variables, will be solved by objective bridges
                "solve_time", "raw_status_string",
                "solve_singlevariable_obj",
                # Quadratic functions are not supported
                "solve_qcp_edge_cases", "solve_qp_edge_cases",
                # Integer and ZeroOne sets are not supported
                "solve_integer_edge_cases", "solve_objbound_edge_cases",
                "solve_zero_one_with_bounds_1",
                "solve_zero_one_with_bounds_2",
                "solve_zero_one_with_bounds_3",
                # Underflow results when using Float64
                "solve_affine_equalto"]
            MOIT.unittest(bridged, config, exclusion_list)
        end
        @testset "Linear tests" begin
            # See explanation in `MOI/test/Bridges/lazy_bridge_optimizer.jl`.
            # This is to avoid `Variable.VectorizeBridge` which does not support
            # `ConstraintSet` modification.
            MOIB.remove_bridge(bridged, MOIB.Constraint.ScalarSlackBridge{T})
            MOIT.contlineartest(bridged, config, [
                # `MOI.UNKNOWN_RESULT_STATUS` instead of `MOI.INFEASIBILITY_CERTIFICATE`
                "linear8a",
                "linear12"
            ])
        end
        @testset "Conic tests" begin
            MOIT.contconictest(bridged, config, [
                # `MOI.UNKNOWN_RESULT_STATUS` instead of `MOI.INFEASIBILITY_CERTIFICATE`
                "lin3", "soc3", "norminf2", "normone2",
                # Missing bridges
                "rootdets",
                # Does not support power and exponential cone
                "pow", "logdet", "exp"])
        end
    end
end
