# Possible Issues and Troubleshooting

We now demonstrate some current limitations of this package via
[`Convex.jl#MathOptInterface`](https://github.com/ericphanson/Convex.jl/tree/MathOptInterface)'s
`Problem Depot`. This is run with `TEST=true`, meaning the solution returned by
the solver will be tested against the true solution.

## Underflows

This occurs when the precision used to represent the solution is not high enough
compared to the internal precision used by the solver. This lack of precision
can lead to catastrophic cancellation. In the following, SDPA-QD is used to
solve the problem, and `Float64` numbers are used to represent the obtained
solution, and the test fails.

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

We try to automatically detect underflows and warn against them; in this case,
the warning is issued. Sometimes this can be avoided by choosing a better set of
parameters. See [Choice of parameters](@ref).

## Presolve

The usefulness of our `presolve` routine is demonstrated in [Using
presolve](@ref). However, our `presolve` subroutine simply uses naïve Gaussian
elimination and has its limitations. At its core, the `reduce!` method takes in
a sparse matrix where each row is a linearized constraint matrix and apply
Gaussian elimination with pivoting to identify the linear independence. This
process is not numerically stable, and we cannot guarantee that all linear
dependent entries will be identified.

This is demonstrated in the following example. We explicitly construct a matrix
with linearly dependent rows. However, due to its numerical instability, we can
only identify a subset of them.

```@repl convex
M = sprand(100, 20000, 0.03);
λ = rand();
M[1, :] = λ*M[3, :] + (1-λ)*M[7, :];
rows = Set(rowvals(SDPAFamily.reduce!(M)[:, 1:end-1]));
redundant = collect(setdiff!(Set(1:10), rows));
@test length(redundant) >= 1
```

## Choice of parameters

Unfortunately, we have not been able to successfully solve every problem that we
have tried with one choice of parameters. We have chosen default parameter
settings that we hope will work with a wide variety of problems. See
[Usage](@ref) for details on switching to two other sets of parameters provided
by the solvers.

This is an example where a better choice of parameters can help.

```@repl convex
test_problem = Convex.ProblemDepot.PROBLEMS["lp"]["lp_dotsort_atom"];
TEST = true; atol = 1e-3; rtol = 0.0;
test_problem(Val(TEST), atol, rtol, Float64) do problem
    solve!(problem, SDPAFamily.Optimizer{Float64}(variant=:sdpa_dd, silent=true))
end
```

```@repl convex
test_problem = Convex.ProblemDepot.PROBLEMS["lp"]["lp_dotsort_atom"];
TEST = true; atol = 1e-3; rtol = 0.0;
test_problem(Val(TEST), atol, rtol, Float64) do problem
    solve!(problem, SDPAFamily.Optimizer{Float64}(variant=:sdpa_dd, silent=true, params = SDPAFamily.UNSTABLE_BUT_FAST))
end
```

## Summary of problematic problems

Due to the above reasons, we have modified the default settings for the
following tests from `Convex.jl`'s `Problem Depot'.

| Solver      | Underflow                                         | Need to use `params = SDPAFamily.UNSTABLE_BUT_FAST`                          | Presolve disabled due to long runtime                  |
| :---------- | :------------------------------------------------ | :----------------------------------------------------------- | :----------------------------------------------------- |
| `:sdpa_dd`  | `affine_Partial_transpose`                        | `affine_Partial_transpose` `lp_pos_atom` `lp_neg_atom` `sdp_matrix_frac_atom` `lp_dotsort_atom` | `affine_Partial_transpose` `lp_min_atom` `lp_max_atom` |
| `:sdpa_qd`  | `affine_Partial_transpose` `affine_Diagonal_atom` | `affine_Partial_transpose` `affine_Diagonal_atom`            | `affine_Partial_transpose` `lp_min_atom` `lp_max_atom` |
| `:sdpa_gmp` | `affine_Partial_transpose`                        | `affine_Partial_transpose`                                   | `affine_Partial_transpose` `lp_min_atom` `lp_max_atom` |

In addition, we have excluded `lp_dotsort_atom` and `lp_pos_atom` when testing
`:sdpa` due to imprecise solutions using default parameters. We have also
excluded all second-order cone problems when using `BigFloat` or `Double64`
numeric types, due to
[MathOptInterface.jl#876](https://github.com/JuliaOpt/MathOptInterface.jl/issues/876),
as well as the `sdp_lambda_max_atom` problem due to
[GenericLinearAlgebra#47](https://github.com/JuliaLinearAlgebra/GenericLinearAlgebra.jl/issues/47).
Both these issues have been fixed on the master branches of the respective
packages, so these exclusions will be removed once new versions are released.

## Troubleshooting

When the solvers fail to return a solution, we recommend trying out the
following troubleshoot steps.

1. Set `silent=false` and look for warnings and error messages. If necessary,
   check the output file. Its path is printed by the solver output and can also
   be retrieved via `Optimizer.tempdir`.
2. Set `presolve=true` to remove redundant constraints. Typically, redundant
   constraints are indicated by a premature `cholesky miss` error as shown
   above.
3. Use `BigFloat` (the default) or `Double64` (from the
   [DoubleFloats](https://github.com/JuliaMath/DoubleFloats.jl) package)
   precision instead of `Float64` (e.g. `SDPAFamily.Optimizer{Double64}(...)`).
   This will reduce the chance of having underflow errors when reading back the
   results.
4. Change the parameters by passing a custom parameter file (i.e.
   `SDPAFamily.Optimizer(params=...)`). [SDPA users
   manual](https://sourceforge.net/projects/sdpa/files/sdpa/sdpa.7.1.1.manual.20080618.pdf)
   contains two other sets of parameters, `UNSTABLE_BUT_FAST` and
   `STABLE_BUT_SLOW`, which can be set by the `params` argument. It might also
   be helpful to use a tighter `epsilonDash` and `epsilonStar` tolerance in a
   custom params file.
