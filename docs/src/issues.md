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

The presolve described in [Presolve](@ref) takes the data as exact numbers of
type `T` and decides with tolerances. The following limitations remain.

* Constraints that are nearly dependent, but not dependent up to the data
  tolerance, are kept. The elimination then divides by small numbers. If the
  mapped-back solution fails the checks, the presolve falls back to the
  unreduced problem.
* Data rounded to a coarser type than `T` (such as `Float64` decimals in a
  `BigFloat` model) satisfy linear relations only up to that rounding; the
  presolve treats such rows as independent unless `presolve_tol` is set (see
  [Arithmetic and tolerances](@ref)).
* A computed number that is genuinely nonzero but smaller than about `τ` times
  the magnitude of its computation (a dynamic range above about `1/τ`, roughly
  `1e31` for `Float64` and `1e125` for `BigFloat` at 256 bits) is treated as zero.
* The scaling equilibrates the constraint matrix and multiplies the right-hand
  side and the objective each by one factor. If the right-hand sides of the
  equilibrated rows differ by many orders of magnitude, for example an inactive
  bound `x ≤ 1e150` next to bounds near 1, one factor cannot bring all of them
  near 1, and the problem can remain too badly scaled for SDPA; the checks then
  report `NUMERICAL_ERROR` instead of a wrong status.
* The checks of the solution returned by SDPA are gross. A solution that is
  close to optimal for a slightly perturbed problem passes them, also when it is
  far from the solution of the given problem because the problem is
  ill-conditioned.
* With `T = Float64`, the working type `Float64x3` has the exponent range of
  `Float64`; for data below about `1e-290` its lower parts are subnormal and the
  rounding tolerance does not hold.
* The fallback drops all reductions, even when only one step of the presolve
  failed its check.
* The eliminations store rows as dictionaries and back-substitute to reduced row
  echelon form. This is fast for sparse equality constraints; for hundreds of
  dense equality constraints with `BigFloat` it takes seconds.

In these cases, `presolve = false` passes the problem to SDPA unchanged.

Two properties of the SDPA interface are independent of the presolve:

* The SDPA status `pdFEAS` (both iterates feasible, iteration limit reached) is
  reported as `OPTIMAL`, even when the duality gap is not small.
* With `presolve = false`, an `INFEASIBILITY_CERTIFICATE` computed by SDPA is
  SDPA's last iterate, which need not be a ray (`A^T y = 0`). With the presolve,
  certificates are checked as described in [Checks and fallback](@ref).

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
