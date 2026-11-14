# Possible Issues and Troubleshooting

We now demonstrate some current limitations of this package via
[`Convex.jl#master`](https://github.com/JuliaOpt/Convex.jl)'s
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
Pkg.add(PackageSpec(name="Convex", url="https://github.com/JuliaOpt/Convex.jl", rev="master"));
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

## Limitations of the presolve

The presolve described in [Presolve](@ref) decides linear dependence with a
tolerance: a row or column counts as dependent if, after Gaussian elimination,
its largest remaining entry is at most `sqrt(eps(T))` times its largest
original entry. Two cases follow from this.

* Rows or columns that are nearly, but not exactly, linearly dependent (for
  example because the problem data were rounded to `Float64` before they were
  converted to `BigFloat`) are kept. SDPA then receives a problem whose Schur
  complement matrix is close to singular.
* The elimination can fill in sparse rows, and the presolve is not optimized for
  speed. In one test with 1000 random sparse equality constraints on a 60×60
  matrix variable, the presolve and writing the input file took about 20 seconds
  with `BigFloat`.
  `presolve = false` skips it.

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

| Solver      | Underflow                                         | Need to use `params = SDPAFamily.UNSTABLE_BUT_FAST`                          |
| :---------- | :------------------------------------------------ | :----------------------------------------------------------- |
| `:sdpa_dd`  | `affine_Partial_transpose`                        | `affine_Partial_transpose` `lp_pos_atom` `lp_neg_atom` `sdp_matrix_frac_atom` `lp_dotsort_atom` |
| `:sdpa_qd`  | `affine_Partial_transpose` `affine_Diagonal_atom` | `affine_Partial_transpose` `affine_Diagonal_atom`            |
| `:sdpa_gmp` | `affine_Partial_transpose`                        | `affine_Partial_transpose`                                   |

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
2. Check that `presolve=true` (the default) is set. Linearly dependent
   constraints typically show up as `cholesky miss condition` warnings. With
   `verbose = SDPAFamily.VERBOSE`, the optimizer prints how many variables the
   presolve eliminated or fixed.
3. Use `BigFloat` (the default) or `Double64` (from the
   [DoubleFloats](https://github.com/JuliaMath/DoubleFloats.jl) package)
   precision instead of `Float64` (e.g. `SDPAFamily.Optimizer{Double64}(...)`).
   This will reduce the chance of having underflow errors when reading back the
   results.
4. Change the parameters by passing choices of parameters, e.g.
    `SDPAFamily.Optimizer(params=(SDPAFamily.UNSTABLE_BUT_FAST))`, or
   `SDPAFamily.Optimizer(params=(epsilonDash=1e-40,))`. It might
   be helpful to use a tighter choice of the `epsilonDash` and `epsilonStar` parameters.
   See [Usage](@ref) and [Changing parameters & solving at very high precision](@ref)
   for more options and examples.
