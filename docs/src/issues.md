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
  tolerance, are kept. The elimination then divides by small numbers. If what
  the presolve removed does not hold at the mapped-back solution, the presolve
  falls back to the unreduced problem.
* Data rounded to a coarser type than `T` (such as `Float64` decimals in a
  `BigFloat` model) satisfy linear relations only up to that rounding; the
  presolve treats such rows as independent unless `presolve_tol` is set (see
  [Arithmetic and tolerances](@ref)). A coefficient that is a rounding error of
  the data (for example `1e-17` where the model means zero) is a coefficient
  like any other.
* A computed number that is genuinely nonzero but more than ``1.5 p + 32`` bits
  below its largest term (a dynamic range of about `1e33` for `Float64` and
  `1e125` for `BigFloat` at 256 bits within one computation) is treated as zero.
* SDPA's tolerances are absolute and apply to the scaled problem. A quantity
  that is small in the scaled problem is not resolved, also when it is large in
  the units of the model:
  - A variable whose whole feasible range is below the tolerance in its scaled
    unit can violate its bounds by that amount (with `1e-310 x + y == 1`,
    `0 <= x <= 1e300`, `y >= 0`, the unit of `x` is `1e310`).
  - A constraint that is violated by less than the tolerance relative to its
    equilibrated row is satisfied for SDPA; an infeasible problem of this kind
    is reported as `OPTIMAL`. This occurs when variables in very different
    units meet in one constraint and the tolerance is `1e-7` (the variant
    `:sdpa`, or `T = Float64`).
  - Ruiz equilibration makes the largest entries of all rows and variables 1.
    For a sparse matrix this does not determine the scaling, and a part of the
    problem can end up at a size below the tolerance of the variant `:sdpa`.
  - If the optimal objective value is far smaller than the product of the
    largest entries of the equilibrated right-hand side and objective, the gap
    tolerance does not resolve it. This case is detected (see
    [Checks and fallback](@ref)) and gives `NUMERICAL_ERROR`; an example is an
    inactive variable with cost 1 that may reach `1e150` on the feasible set.
* SDPA reports unboundedness when the objective of the scaled problem passes
  its parameters `lowerBound` or `upperBound` (`±1e5` by default). A problem
  whose optimal value in the scaled problem is beyond these bounds is reported
  as `DUAL_INFEASIBLE` or `INFEASIBLE`.
* The checks of the solution returned by SDPA are at the square root of SDPA's
  tolerances. A solution that is close to optimal for a slightly perturbed
  problem passes them, also when it is far from the solution of the given
  problem because the problem is ill-conditioned.
* The fallback drops all reductions, even when only one step of the presolve
  failed its check.
* The equality elimination computes the reduced constraints explicitly. For
  `m` dense equality constraints this takes time proportional to `m^3` in
  `BigFloat` arithmetic (about 6 seconds for `m = 400` at 256 bits).

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
