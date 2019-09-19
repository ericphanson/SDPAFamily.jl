using Convex, SDPAFamily, Test, MathOptInterface
const MOI = MathOptInterface


problem_func = Convex.ProblemDepot.PROBLEMS["lp"]["lp_dotsort_atom"];
problem_func(Val(true), 1e-3, 0.0, BigFloat) do p
    opt = SDPAFamily.Optimizer{}(; 
    params = SDPAFamily.UNSTABLE_BUT_FAST,
    presolve = false, variant = :sdpa_dd)
    Convex.solve!(p, opt)
    @test MOI.get(opt, MOI.TerminationStatus()) == MOI.ITERATION_LIMIT
end
