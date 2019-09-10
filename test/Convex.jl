using Convex: Convex
using Convex.ProblemDepot: run_tests, foreach_problem


@testset "Convex tests with variant $var" begin
    foreach_problem(;exclude=[  r"mip",
                                r"exp",
                                r"benchmark",
                                r"lp_min_atom", # hang
                            ]) do name, problem_func
        @testset "$name" begin
            problem_func(Val(true), 1e-3, 0.0, Float64) do p
                @info "Testing" name var
                solve!(p, SDPAFamily.Optimizer{Float64}(presolve = true, silent = true, variant = var))
            end
        end
    end
end
