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
    # These problems use random data. With the presolve, some solves fail the
    # tests for `:sdpa_dd`, `:sdpa_qd` and `:sdpa_gmp`, and the failures change
    # from run to run.
    r"sdp_quantum_relative_entropy1_lowrank",
    r"sdp_quantum_relative_entropy2_lowrank",
    r"sdp_quantum_relative_entropy3_lowrank",
    r"sdp_quantum_relative_entropy4_lowrank",
    # With the presolve, this passes for `:sdpa`, `:sdpa_qd` and `:sdpa_gmp`, but
    # each solve takes about a minute; for `:sdpa_dd` the optimal value is off by
    # more than the tolerance in 5 of 6 runs.
    r"sdp_lieb_ando",
]

# Problems that cannot be handled due to issues with a certain numeric type,
# independent of variant.
const type_excludes = Dict( Float64 => Regex[],
                            BigFloat => Regex[
                                # `DomainError`: "`A` cannot have complex eigenvalues."
                                # when the test computes the eigenvalues of the solution.
                                # See https://github.com/jump-dev/Convex.jl/issues/502
                                r"sdp_dual_lambda_max_atom",
                                r"sdp_lambda_min_atom",
                            ],
                            Double64 => Regex[])

# Problems that cannot be handled with a specific combination of variant and numeric type
# (often due to things like underflow).
const variant_excludes = Dict(
                    (:sdpa, BigFloat) => Regex[
                            r"sdp_trace_logm_real", # failed 1 of 2 runs (imprecise)
                        ],
                    )
for T in (Float64, Double64, BigFloat), var in (:sdpa_qd, :sdpa_gmp)
    # wrong optimal value in some runs
    push!(get!(variant_excludes, (var, T), Regex[]), r"sdp_quantum_channel_capacity")
end

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
                            settings = (presolve = true, variant = var)

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
