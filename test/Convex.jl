using Test, SDPAFamily, DoubleFloats
using Convex: Convex
using Convex.ProblemDepot: run_tests, foreach_problem
using GenericLinearAlgebra

# Problems that cannot be handled by any variant with any numeric type
const common_excludes = Regex[
    r"mip", # SDPA solvers don't support mixed integer programming
    r"exp", # SDPA solvers don't support the exponential cone (no bridge yet?)
    r"benchmark", # don't include benchmark-only problems
    r"sdp_Complex_Semidefinite_constraint", # too large / slow
    # The following were added in https://github.com/ericphanson/SDPAFamily.jl/pull/66
    # They fail because of issues like https://github.com/jump-dev/Convex.jl/issues/502 or
    # the fact that we removed the presolve, see https://github.com/jump-dev/Convex.jl/issues/503
    r"lp_dotsort_atom",
    r"lp_sumlargest_atom",
    r"sdp_Real_Variables_with_complex_equality_constraints",
    r"sdp_dual_lambda_max_atom",
    r"sdp_geom_mean_epicone_real_8_5_optA",
    r"sdp_geom_mean_epicone_real_neg3_5_optA",
    r"sdp_lambda_min_atom",
    r"sdp_lieb_ando",
    r"sdp_quantum_relative_entropy2_lowrank",
    r"sdp_quantum_relative_entropy3_lowrank",
    r"sdp_sdp_constraints",
    r"sdp_sum_largest_eigs",
    r"sdp_trace_logm_cplx",
    r"sdp_trace_logm_real",
    r"sdp_trace_mpower_cplx_2_3",
    r"sdp_trace_mpower_cplx_5_4",
    r"sdp_trace_mpower_cplx_neg1_4",
    r"sdp_trace_mpower_real_2_3",
    r"sdp_trace_mpower_real_5_4",
    r"sdp_trace_mpower_real_neg1_4",
]

# Problems that cannot be handled due to issues with a certain numeric type,
# independent of variant.
const type_excludes = Dict( Float64 => Regex[],
                            BigFloat => Regex[],
                            Double64 => Regex[])

# Problems that cannot be handled with a specific combination of variant and numeric type
# (often due to things like underflow).
const variant_excludes = Dict(
                    (:sdpa, Float64) => Regex[
                            r"lp_dotsort_atom", # imprecise, cholesky miss
                            r"lp_pos_atom", # imprecise
                            r"sdp_matrix_frac_atom", # imprecise
                        ],
                    (:sdpa, BigFloat) => Regex[
                            r"lp_dotsort_atom", # imprecise, cholesky miss
                            r"lp_pos_atom", # imprecise
                            r"sdp_matrix_frac_atom", # imprecise
                        ],
                    (:sdpa, Double64) => Regex[
                            r"lp_dotsort_atom", # imprecise, cholesky miss
                            r"lp_pos_atom", # imprecise
                            r"sdp_matrix_frac_atom", # imprecise
                        ]
                    )

# failures on CI with Julia 1.0; couldn't locally reproduce locally
if VERSION < v"1.3"
    for T in (Float64, Double64, BigFloat)
        for variant in variants
            variant_excludes[(variant, T)] = Regex[r"socp_quad_form_atom"]
        end
    end

    for T in (Float64, Double64, BigFloat)
        variant_excludes[(:sdpa, T)] = Regex[r"lp_pos_atom"]
    end
end

# problems where `presolve=true` causes problems
const no_presolve_problems = ["affine_Partial_transpose", "lp_min_atom", "lp_max_atom"]


# Some problems need a different choice of parameters to pass the tests
const params_options = Dict(
                    (:sdpa_gmp, Float64) => Dict(
                            "affine_Partial_transpose" => SDPAFamily.UNSTABLE_BUT_FAST,
                        ),
                    (:sdpa_gmp, Double64) => Dict(
                            "affine_Partial_transpose" => SDPAFamily.UNSTABLE_BUT_FAST,
                        ),
                    (:sdpa_dd, Float64) =>  Dict(
                            "lp_dotsort_atom" => SDPAFamily.UNSTABLE_BUT_FAST,
                            "lp_pos_atom" => SDPAFamily.UNSTABLE_BUT_FAST,
                            "lp_neg_atom" => SDPAFamily.UNSTABLE_BUT_FAST,
                            "sdp_matrix_frac_atom" => SDPAFamily.UNSTABLE_BUT_FAST,
                            "affine_Partial_transpose" => SDPAFamily.UNSTABLE_BUT_FAST,
                            "affine_Diagonal_atom" => SDPAFamily.UNSTABLE_BUT_FAST,
                        ),
                    (:sdpa_dd, BigFloat) =>  Dict(
                        "lp_dotsort_atom" => SDPAFamily.UNSTABLE_BUT_FAST,
                        "lp_pos_atom" => SDPAFamily.UNSTABLE_BUT_FAST,
                        "lp_neg_atom" => SDPAFamily.UNSTABLE_BUT_FAST,
                        "sdp_matrix_frac_atom" => SDPAFamily.UNSTABLE_BUT_FAST,
                        "affine_Partial_transpose" => SDPAFamily.UNSTABLE_BUT_FAST,
                        "affine_Diagonal_atom" => SDPAFamily.UNSTABLE_BUT_FAST,
                        ),
                    (:sdpa_dd, Double64) =>  Dict(
                            "lp_dotsort_atom" => SDPAFamily.UNSTABLE_BUT_FAST,
                            "lp_pos_atom" => SDPAFamily.UNSTABLE_BUT_FAST,
                            "lp_neg_atom" => SDPAFamily.UNSTABLE_BUT_FAST,
                            "sdp_matrix_frac_atom" => SDPAFamily.UNSTABLE_BUT_FAST,
                            "affine_Partial_transpose" => SDPAFamily.UNSTABLE_BUT_FAST,
                            "affine_Diagonal_atom" => SDPAFamily.UNSTABLE_BUT_FAST,
                            ),
                    (:sdpa_qd, Float64) =>  Dict(
                            "affine_Partial_transpose" => SDPAFamily.UNSTABLE_BUT_FAST,
                            "affine_Diagonal_atom" => SDPAFamily.UNSTABLE_BUT_FAST,
                        ),
                    (:sdpa_qd, BigFloat) =>  Dict(
                            "affine_Partial_transpose" => SDPAFamily.UNSTABLE_BUT_FAST,
                            "affine_Diagonal_atom" => SDPAFamily.UNSTABLE_BUT_FAST,
                        ),
                    (:sdpa_qd, Double64) =>  Dict(
                            "affine_Partial_transpose" => SDPAFamily.UNSTABLE_BUT_FAST,
                            "affine_Diagonal_atom" => SDPAFamily.UNSTABLE_BUT_FAST,
                        ))
@testset "Convex tests" for var in variants
    @testset "Convex tests with variant $var and type $T" for T in (Float64, Double64, BigFloat)
    @info "Starting testset `Convex tests with variant $var and type $T`"

        excludes = vcat(common_excludes, get(variant_excludes, (var, T), Regex[]), type_excludes[T])
        for class in keys(Convex.ProblemDepot.PROBLEMS)
            @testset "$class" begin
                foreach_problem(class; exclude = excludes) do name, problem_func
                    @testset "$name" begin
                        problem_func(Val(true), 1e-3, 0.0, T) do p
                            # @info "`solve!` called" name var T
                            presolve = !(name ∈ no_presolve_problems)
                            settings = (presolve = presolve, variant = var)

                            params = get(get(params_options, (var, T), Dict()), name, nothing)
                            if params !== nothing
                                settings = (params = params, settings...)
                            end
                            @info "Solving problem $name with variant $var and type $T"
                            Convex.solve!(p, () -> SDPAFamily.Optimizer{T}(; settings...))
                        end
                    end
                end
            end
        end
    end
end
