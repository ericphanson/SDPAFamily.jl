using Convex: Convex
using Convex.ProblemDepot: run_tests, foreach_problem

const common_excludes = [
                            r"mip", # doesn't support MIP
                            r"exp", # doesn't support exponential cone
                            r"benchmark", # don't include benchmark-only problems
                            r"sdp_Complex_Semidefinite_constraint", # too large
                        ]

const excludes_dict = Dict(:sdpa_gmp => vcat(common_excludes,
                        Regex[r"affine_Partial_transpose", # underflows

                        ]),
                    :sdpa_dd =>  vcat(common_excludes,
                        Regex[r"affine_Partial_transpose", # underflows

                        ]),
                    :sdpa_qd =>  vcat(common_excludes,
                        Regex[
                            r"affine_Partial_transpose", # underflows
                            r"affine_Diagonal_atom" # underflows
                        ]),
                    :sdpa => vcat(common_excludes,
                        Regex[
                            r"lp_dotsort_atom", # imprecise
                            r"lp_pos_atom" # imprecise
                        ]))

const no_presolve_problems = ["affine_Diagonal_atom", "affine_Partial_transpose", "lp_min_atom", "lp_max_atom"]

@testset "Convex tests with variant $var" begin
    foreach_problem(;exclude=excludes_dict[var]) do name, problem_func
        @testset "$name" begin
            problem_func(Val(true), 1e-3, 0.0, Float64) do p
                @info "`solve!` called" name var
                presolve = !(name ∈ no_presolve_problems)
                time = @elapsed solve!(p, SDPAFamily.Optimizer{Float64}(presolve = presolve, silent = true, variant = var))
                @info "Finished `solve!`" time
            end
        end
    end
end
