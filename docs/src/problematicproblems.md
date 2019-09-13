# Problematic Problems

Unfortunately, we have not been able to successfully solve every problem that we have tried with one choice of parameters. We have chosen default parameter settings that we hope will work with a wide variety of problems.

## Underflows

This occurs when the precision used to represent the solution is not high enough compared to the internal precision used by the solver. This lack of precision can lead to catastraphic cancellation. For example, the following problem is taken from the Convex.jl Problem Depot, and is run with `TEST=true`, meaning the solution returned by the solver will be tested against the true solution. In the following, SDPA-QD is used to solve the problem, and `Float64` numbers are used to represent the obtained solution, and the test fails.

```@repl
using SDPAFamily
using Pkg
Pkg.add(PackageSpec(name="Convex", url="https://github.com/ericphanson/Convex.jl", rev="MathOptInterface"));
using Convex

test_problem = Convex.ProblemDepot.PROBLEMS["affine"]["affine_Diagonal_atom"]
TEST = true; atol = 1e-3; rtol = 0.0
test_problem(Val(TEST), atol, rtol, Float64) do problem
    solve!(problem, SDPAFamily.Optimizer{Float64}(variant=:sdpa_qd, silent=true))
end
```

We try to automatically detect underflows and warn against them; in this case, the warning is issued.

## Presolve

[todo: add example where presolve helps and an example where presolve backfires]

## Troubleshooting

Try:

1. Set `silent=false` and look for warnings and error messages
2. Set `presolve=true` to remove redundant constraints
3. Use `BigFloat` (the default) or `Double64` (from the [DoubleFloats](https://github.com/JuliaMath/DoubleFloats.jl) package) precision instead of `Float64` (e.g. `SDPAFamily.Optimizer{Double64}(...)`)
4. Change the parameters by passing a custom parameter file (i.e. `SDPAFamily(params_path=...)`).
