using Test, SDPAFamily
using Convex: Convex
using Convex.ProblemDepot: run_tests, foreach_problem
using GenericLinearAlgebra

# Problems that cannot be handled by any variant with any numeric type
const common_excludes = Regex[
                            r"mip", # SDPA solvers don't support mixed integer programming
                            r"exp", # SDPA solvers don't support the exponential cone (no bridge yet?)
                            r"benchmark", # don't include benchmark-only problems
                            r"sdp_Complex_Semidefinite_constraint", # too large / slow
                        ]

# Problems that cannot be handled due to issues with a certain numeric type,
# independent of variant.
const type_excludes = Dict( Float64 => Regex[],
                            BigFloat => Regex[
                                r"sdp_lambda_max_atom", # GenericLinearAlgebra#47
                                r"socp", # MathOptInterface.jl#876
                        ])

# Problems that cannot be handled with a specific combination of variant and numeric type
# (often due to things like underflow).
const variant_excludes = Dict(
                    (:sdpa_dd, BigFloat) =>  Regex[
                            r"affine_Diagonal_atom", # needs smaller epsilon
                        ],
                    (:sdpa_qd, BigFloat) =>  Regex[
                            r"affine_Diagonal_atom", # needs smaller epsilon
                        ],
                    (:sdpa, Float64) => Regex[
                            r"lp_dotsort_atom", # imprecise, cholesky miss
                            r"lp_pos_atom" # imprecise
                        ],
                    (:sdpa, BigFloat) => Regex[
                            r"lp_dotsort_atom", # imprecise, cholesky miss
                            r"lp_pos_atom" # imprecise
                        ])

# problems where `presolve=true` causes problems
const no_presolve_problems = ["affine_Partial_transpose", "lp_min_atom", "lp_max_atom"]


# Some problems need a different choice of parameters to pass the tests
const params_options = Dict(
                    (:sdpa_gmp, Float64) => Dict(
                            "affine_Partial_transpose" => "-pt 1",
                        ),
                    (:sdpa_dd, Float64) =>  Dict(
                            "lp_dotsort_atom" => "-pt 1",
                            "lp_pos_atom" => "-pt 1",
                            "lp_neg_atom" => "-pt 1",
                            "sdp_matrix_frac_atom" => "-pt 1",
                            "affine_Partial_transpose" => "-pt 1",
                        ),
                    (:sdpa_dd, BigFloat) =>  Dict(
                        "lp_dotsort_atom" => "-pt 1",
                        "lp_pos_atom" => "-pt 1",
                        "lp_neg_atom" => "-pt 1",
                        "sdp_matrix_frac_atom" => "-pt 1",
                        "affine_Partial_transpose" => "-pt 1",
                        ),
                    (:sdpa_qd, Float64) =>  Dict(
                            "affine_Partial_transpose" => "-pt 1",
                            "affine_Diagonal_atom" => "-pt 1",
                        ),
                    (:sdpa_qd, BigFloat) =>  Dict(
                            "affine_Partial_transpose" => "-pt 1",
                        ))

@testset "Convex tests with variant $var and type $T" for T in (Float64, BigFloat)
    excludes = vcat(common_excludes, get(variant_excludes, (var, T), Regex[]), type_excludes[T])
    foreach_problem(;exclude = excludes) do name, problem_func
        @testset "$name" begin
            problem_func(Val(true), 1e-3, 0.0, T) do p
                # @info "`solve!` called" name var T
                presolve = !(name ∈ no_presolve_problems)
                settings = (presolve = presolve, silent = true, variant = var)

                params = get(get(params_options, (var, T), Dict()), name, nothing)
                if params !== nothing
                    settings = (params_path = params, settings...)
                end

                Convex.solve!(p, SDPAFamily.Optimizer{T}(; settings...))

            end
        end
    end
end
