# Problematic Problems

Unfortunately, we have not been able to successfully solve every problem that we have tried with one choice of parameters. We have chosen default parameter settings that we hope will work with a wide variety of problems.

We now demonstrate some common problems using [`Convex.jl#MathOptInterface`](https://github.com/ericphanson/Convex.jl/tree/MathOptInterface).

## Underflows

This occurs when the precision used to represent the solution is not high enough compared to the internal precision used by the solver. This lack of precision can lead to catastrophic cancellation. For example, the following problem is taken from the `Convex.jl Problem Depot`, and is run with `TEST=true`, meaning the solution returned by the solver will be tested against the true solution. In the following, SDPA-QD is used to solve the problem, and `Float64` numbers are used to represent the obtained solution, and the test fails.

```@setup convex
using SDPAFamily, Test, SparseArrays
using Pkg
Pkg.add(PackageSpec(name="Convex", url="https://github.com/ericphanson/Convex.jl", rev="MathOptInterface"));
using Convex
```
```@repl convex
test_problem = Convex.ProblemDepot.PROBLEMS["affine"]["affine_Diagonal_atom"];
TEST = true; atol = 1e-3; rtol = 0.0;
test_problem(Val(TEST), atol, rtol, Float64) do problem
    solve!(problem, SDPAFamily.Optimizer{Float64}(variant=:sdpa_qd, silent=true))
end
```

We try to automatically detect underflows and warn against them; in this case, the warning is issued.

## Presolve

The usefulness of our `presolve` routine is demonstrated in [Using presolve](@ref). However, our `presolve` subroutine does not utilize any tools more powerful than naïve Gaussian elimination and it has limitations. At its core, the `reduce!` method takes in a sparse matrix where each row is a linearized constraint matrix and apply Gaussian elimination with pivoting to identify the linear independence. This process is not numerically stable as the rounding errors accumulate as we progress. Therefore, we cannot guarantee that all linear dependent entries will be identified. 

This is demonstrated in the following example. For a ``1500 \times 100`` matrix (Here we use one more column because `reduce!` ignores the last column as it represents the constraint constants. They go along for the ride here to help us identify _inconsistent_ constraints.), we expect ``1400`` or more linearly dependent rows. However, due to its numerical instability, we can only identify a subset of them. 

```@repl convex
A = sprandn(1500, 101, 0.1);
redundant_F = collect(setdiff!(Set(1:150), Set(rowvals(SDPAFamily.reduce!(A)[:, 1:end-1]))));
length(redundant_F)
```

## Troubleshooting

Try:

1. Set `silent=false` and look for warnings and error messages
2. Set `presolve=true` to remove redundant constraints
3. Use `BigFloat` (the default) or `Double64` (from the [DoubleFloats](https://github.com/JuliaMath/DoubleFloats.jl) package) precision instead of `Float64` (e.g. `SDPAFamily.Optimizer{Double64}(...)`)
4. Change the parameters by passing a custom parameter file (i.e. `SDPAFamily(params_path=...)`).
