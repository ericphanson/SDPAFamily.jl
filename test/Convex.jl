using Convex: Convex
using Convex.ProblemDepot: run_tests, foreach_problem

common_excludes = [  r"mip",
                    r"exp",
                    r"benchmark",
                    r"lp_min_atom", # hangs
                    r"lp_max_atom", # hangs
                    ]
excludes_dict = Dict(:gmp => vcat(common_excludes,Regex[
                        r"affine_Partial_transpose" # slow
                    ]),
                    :dd =>  vcat(common_excludes,Regex[
                        r"affine_Partial_transpose" # slow
                    ]),
                    :qd =>  vcat(common_excludes,Regex[
                        r"affine_Partial_transpose" # slow
                    ]),
                    :plain => vcat(common_excludes, Regex[]),)

@testset "Convex tests with variant $var" begin
    foreach_problem(;exclude=excludes_dict[var]) do name, problem_func
        @testset "$name" begin
            problem_func(Val(true), 1e-3, 0.0, Float64) do p
                @info "`solve!` called" name var
                time = @elapsed solve!(p, SDPAFamily.Optimizer{Float64}(presolve = true, silent = true, variant = var))
                @info "Finished `solve!`" time
            end
        end
    end
end
