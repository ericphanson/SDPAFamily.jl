using Convex: Convex
using Convex.ProblemDepot: run_tests


@testset "Convex tests with variant $var" begin
    run_tests(; exclude=[   r"mip",
                            r"socp"
                        ]) do p
        solve!(p, SDPAFamily.Optimizer{Float64}(presolve = true, silent = true, variant = var))
    end
end
