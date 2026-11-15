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

The [Presolve](@ref) checks every result against the model and reports
`NUMERICAL_ERROR` when a check fails. The cases below are where that happens, or
where a result can be misleading.

* If the variables or constraints of a model are in very different units, the
  result can be `NUMERICAL_ERROR`. SDPA's tolerances are absolute. The presolve
  rescales the problem so that coefficients are near 1, but it cannot always
  tell which of two very different magnitudes in one constraint matters. In
  tests with random problems whose variables were scaled by factors between
  `1e-6` and `1e6` and whose constraints were scaled by factors between `1e-12`
  and `1e12`, about 2% of the feasible problems ended this way with
  `variant = :sdpa` (tolerance `1e-7`) and none with `:sdpa_dd`, `:sdpa_qd` and
  `:sdpa_gmp`. With factors of `1e±30` a few statuses were wrong for all
  variants. Choose units so that variables and constraints have comparable
  sizes.
* A variable that is small in the units of the model, such as `0 <= x <= 1e-9`,
  is solved relative to its bounds if the presolve finds them. It finds bounds
  written as constraints on one variable, and bounds that follow from a linear
  inequality and the bounds of its other variables (`x1 + x2 <= 1e-9` with
  `x1 >= 0` and `x2 >= 0`). Otherwise values are accurate to SDPA's absolute
  tolerance only: the status can be `NUMERICAL_ERROR`, or a result reported as
  `OPTIMAL` can violate a small bound by a large factor of that bound. This
  remains the case for a small range that follows only from a matrix constraint,
  for variables that inequalities force to zero (`x >= 0` together with
  `x <= 0`, directly or through other constraints) and variables whose bounds
  depend on them, for a variable whose bounds
  differ by many orders of magnitude, such as `1e-9 <= x <= 1e9` (the smaller
  bound is resolved relative to the larger one), and for entries of very
  different sizes in one matrix constraint, which is rescaled by a single
  factor. `presolve = false` does not rescale at all.
* If the constraints of a model contradict each other by less than SDPA's
  tolerance, relative to the sizes of the terms in those constraints, the model
  is reported as `OPTIMAL`.
* If the optimal value is very large compared to the data, a feasible and
  bounded model can be reported as `INFEASIBLE` or `DUAL_INFEASIBLE`. SDPA
  reports a problem as unbounded when an objective value passes its parameters
  `lowerBound` or `upperBound`, which are `-1e5` and `1e5` by default. The
  presolve rescales the objective, so this depends on the rescaled problem. An
  example is the chain `x[1] >= 1`, `x[i+1] >= 2 x[i]` for 20 variables with
  objective `sum(x)`, whose optimal value is about `1e6`: with `:sdpa_dd` it is
  solved with the presolve and reported as `INFEASIBLE_OR_UNBOUNDED` with `presolve = false`.
  If a model is affected, set larger bounds, for example
  `params = (lowerBound = -1e10, upperBound = 1e10)`.
* A constraint that follows from others only approximately, with an error above
  `presolve_tol`, is kept. The reduced problem is then ill-conditioned, and the
  result can fail the checks. If the data were rounded to a coarser type than
  `T`, set `presolve_tol` (see [Options](@ref)).
* A coefficient that is a rounding error of the data, such as `1e-17` where the
  model means zero, is a coefficient like any other.
* If a value of the solution is larger than the largest number of type `T`
  (about `1.8e308` for `Float64`), the status is `NUMERICAL_ERROR`.
* The time of the presolve grows with the cube of the number of dense equality
  constraints (see [Cost](@ref)).

With `presolve = false`, the model is passed to SDPA unchanged. An
`INFEASIBILITY_CERTIFICATE` is then SDPA's last iterate and is not checked.

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
