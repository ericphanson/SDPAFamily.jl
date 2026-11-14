# Usage

The main object of interest supplied by this package is
`SDPAFamily.Optimizer{T}()`. Here, `T` is a numeric type which defaults to
`BigFloat`. Keyword arguments may also be passed to `SDPAFamily.Optimizer()`:

* `variant`: either `:sdpa_gmp`, `:sdpa_qd`, or `:sdpa_dd`, to use SDPA-GMP,
  SDPA-QD, or SDPA-DD, respectively. Defaults to `:sdpa_gmp`.
* `silent`: a boolean to indicate whether or not to print output. Defaults to
  `false`.
* `verbose`: accepts either `SDPAFamily.SILENT` (which is equivalent to
  `silent=true`), `SDPAFamily.WARN` (which is the default), or
  `SDPAFamily.VERBOSE`, which prints all output from the binary as well as
  additional warnings.
* `presolve`: whether to run the presolve described in [Presolve](@ref) before
  calling the solver. Defaults to `true`.
* `presolve_tol`: how exactly a constraint must follow from others to count as
  redundant, see [Options](@ref). Defaults to `nothing`, which means
  `64 * eps(T)`.
* `binary_path`: a string representing a path to the SDPA-GMP binary to use. The
  default is chosen at `build` time. To change the default binaries, see [Custom
  binary](@ref).
* `params`: either `SDPAFamily.DEFAULT` (the default option),
  `SDPAFamily.UNSTABLE_BUT_FAST`, `SDPAFamily.STABLE_BUT_SLOW`, a `NamedTuple`
  giving a list of choices of parameters (e.g. `params = (maxIteration=600,)`),
  an [`SDPAFamily.Params`](@ref) object, or or a string representing a path to a
  custom parameter file. See [`SDPAFamily.Params`](@ref) for the possible choices of parameters.

The default parameters used by `SDPAFamily.jl` depend on the variant and numeric
type, and can be found as the default values in the source code
[here](https://github.com/ericphanson/SDPAFamily.jl/tree/master/src/params.jl).
The choices choices of parameters `UNSTABLE_BUT_FAST` and `STABLE_BUT_SLOW` are
documented in the [SDPA users
manual](https://sourceforge.net/projects/sdpa/files/sdpa/sdpa.7.1.1.manual.20080618.pdf),
as well as the meaning of each parameter. A new parameter file is generated for
each optimizer instance `opt` (and is regenerated for each solve), and can be
found in the directory `opt.tempdir`, along with the SDPA-format input and
output files for the problem, assuming a custom parameter file has not been
passed to the optimizer. See [Changing parameters & solving at very high
precision](@ref) below for an example.

`SDPAFamily.Optimizer()` also accepts `variant = :sdpa` to use the
non-high-precision SDPA binary. For general usage of the SDPA solver, use
[SDPA.jl](https://github.com/JuliaOpt/SDPA.jl) which uses the C++ library to
interface directly with the SDPA binary.

## Using a number type other than `BigFloat`

`SDPAFamily.Optimizer()` uses `BigFloat` for problem data and solution by
default. To use, for example, `Float64` instead, simply call
`SDPAFamily.Optimizer{Float64}()`. However, this may cause underflow errors when
reading the solution file, particularly when `MathOptInterface` bridges are
used. Note that with `MathOptInterface` all the problem data must be
parametrised by the same number type, i.e. `Float64` in this case.

## Presolve

SDPA's algorithm needs a problem without equality constraints and without
linearly dependent constraints. Models written with JuMP or Convex.jl usually
have both: equality constraints such as `tr(X) == 1`, and redundant constraints
that the modelling layer adds when it rewrites the model. On such models SDPA
stops at its iteration limit, prints `cholesky miss condition` warnings, or
returns an inaccurate point.

The presolve removes these constraints before SDPA is called. It uses the
equality constraints to eliminate variables, drops constraints that follow from
others, solves the smaller problem with SDPA, and computes the solution of the
model from the result. It is on by default.

```@repl presolve
using SDPAFamily, Convex
x = Variable(); y = Variable();
problem = minimize(x + y, [x + y >= 1]; numeric_type = BigFloat);
solve!(problem, () -> SDPAFamily.Optimizer(presolve = false, silent = true))
problem.status
solve!(problem, () -> SDPAFamily.Optimizer(silent = true))
problem.status
problem.optval
```

In this problem the two variables have the same coefficients. Without the
presolve, SDPA-GMP stops at the iteration limit. With the presolve, one of the
two variables is fixed to zero and SDPA-GMP finds the optimal value 1.

### What you get back

Variable values, constraint values and dual values are those of your model, not
of the smaller problem. A redundant equality constraint that was dropped has the
dual value zero. A variable that was fixed has the value zero.

Some models are decided by the presolve, and SDPA is not called:

| Situation | `TerminationStatus` | Results |
|:--|:--|:--|
| Equality constraints contradict each other | `INFEASIBLE` | The dual solution is a certificate of infeasibility. |
| The objective decreases along a direction that changes no constraint | `DUAL_INFEASIBLE` | The primal solution is an unbounded ray. |
| The equality constraints determine all variables | `OPTIMAL`, or `INFEASIBLE` if another constraint is violated at that point | The point, and dual values if it is optimal. |

`MOI.RawStatusString` says when a result comes from the presolve.

After SDPA has run, the presolve checks the result against your model. Three
outcomes need attention:

* `NUMERICAL_ERROR`: SDPA reported a result, but it does not satisfy the model to
  about half of the digits of SDPA's tolerances, also after a second solve with
  a different rescaling of the problem. The primal and dual statuses are
  `UNKNOWN_RESULT_STATUS`; the values are still available, and
  `MOI.RawStatusString` names SDPA's status and the check that failed. This
  usually means that the model is badly scaled (see
  [Limitations of the presolve](@ref)). Rescale the model, use a variant with
  more precision, or try `presolve = false`.
* `ITERATION_LIMIT`: SDPA stopped at `maxIteration`. If the primal and dual
  statuses are `FEASIBLE_POINT`, both points are feasible but their objective
  values do not agree to half of the digits of SDPA's gap tolerance. Increase
  `maxIteration` or see [Choice of parameters](@ref).
* A warning that the problem "is passed to the solver without the reductions of
  the presolve": the presolve could not confirm its own reduction, for example
  because constraints are almost but not exactly redundant, and SDPA solved the
  model as with `presolve = false`. The result is checked in the same way. If
  the data of the model were rounded, see `presolve_tol` below.

### Options

`presolve = false` passes the model to SDPA unchanged, with each equality
constraint written as two inequalities. Use it to compare, or if the presolve
takes too long on a model that SDPA can solve directly.

`presolve_tol` sets how exactly a constraint must follow from others to count
as redundant. The default is `64 * eps(T)` for the number type `T` of the
optimizer, relative to the sizes of the coefficients involved. The data are
taken as exact numbers of type `T`. This matters when a `BigFloat` model is
built from `Float64` or decimal data: `0.3` and `3 * 0.1` are different
`Float64` numbers, so constraints that are redundant on paper are independent in
the model. In the following problem the second equality constraint is three
times the first on paper. With the default tolerance the two constraints
determine `x`. With `presolve_tol = 64 * eps(Float64)` the second is dropped as
redundant.

```@repl presolve
x = Variable(2);
problem = minimize(x[1] + x[2], [0.1 * x[1] + 0.3 * x[2] == 0.4, 0.3 * x[1] + 0.9 * x[2] == 1.2, x >= 0]; numeric_type = BigFloat);
solve!(problem, () -> SDPAFamily.Optimizer(silent = true))
Float64(problem.optval)
solve!(problem, () -> SDPAFamily.Optimizer(silent = true, presolve_tol = 64 * eps(Float64)))
Float64(problem.optval)
```

`verbose = SDPAFamily.VERBOSE` prints how many variables the presolve eliminated
and fixed, and how long it took.

### Cost

The presolve is Gaussian elimination in extended precision. On an Apple M-series
laptop it took about 0.1 seconds for 1000 sparse equality constraints on a 60×60
matrix variable, and for 2000 sparse equality constraints on 4000 nonnegative
variables. For 400 dense equality constraints on 420 variables it took about 5
seconds. Eliminating variables can make sparse constraints denser: in a test
with 400 variables and 200 sparse equality constraints the SDPA input file was
86 KiB with the presolve and 34 KiB without it.

How the presolve works is described in the comments of
[`src/presolve.jl`](https://github.com/ericphanson/SDPAFamily.jl/blob/master/src/presolve.jl).

## Changing parameters & solving at very high precision

The following example uses the problem from [Optimal guessing probability for a
pair of quantum states](@ref). We can increase the precision by changing the
default parameters:

```@setup convexquantum
using SDPAFamily, Test, Convex, LinearAlgebra
```

```@repl convexquantum
opt = () -> SDPAFamily.Optimizer(
    params = (  epsilonStar = 1e-200, # constraint tolerance
                epsilonDash = 1e-200, # normalized duality gap tolerance
                precision = 2000 # arithmetric precision used in sdpa_gmp
    ))

setprecision(2000) # set Julia's global BigFloat precision to 2000
ρ₁ = Complex{BigFloat}[1 0; 0 0]
ρ₂ = (1//2)*Complex{BigFloat}[1 -im; im 1]
E₁ = ComplexVariable(2, 2);
E₂ = ComplexVariable(2, 2);
problem = maximize( real((1//2)*tr(ρ₁*E₁) + (1//2)*tr(ρ₂*E₂)),
                    [E₁ ⪰ 0, E₂ ⪰ 0, E₁ + E₂ == Diagonal(ones(2))];
                    numeric_type = BigFloat );
solve!(problem, opt)
p_guess = 1//2 + 1/(2*sqrt(big(2)))
@test problem.optval ≈ p_guess atol=1e-200
```

With these parameters, we have recovered the true answer to a tolerance of $10^{-200}$.

Note that we called `setprecision(2000)` at the start. This is so that the `BigFloat` objects used to construct `ρ₁`, `ρ₂`, the internals of the `Convex.Problem` instance are constructed to such a precision, as well as the BigFloat objects used to store the output of SDPA-GMP. The default, `256`, is not sufficient in this case. Testing this can be a little bit subtle: for example, if `setprecision(2000)` was not called (or equivalently, `setprecision(256)` called instead), then the test

```julia
@test problem.optval ≈ p_guess atol=1e-200
```

would still pass. However, that is because `p_guess` is only constructed at approximately 77 digits of precision, and `problem.optval` is only read back from SDPA-GMP at the same precision. So in that case, the test isn't truly testing that the solution is accurate to 200 digits of precision.

In this case, since `1` and `1/2` are exactly representable by floating point numbers, it is enough to specify `setprecision(2000)` before the `solve!` call (so the `ρ₁` and `ρ₂` are only constructed at 256 bits of high precision), but it is good practice to set the precision at the start for the whole problem. Moreover, since the precision is mutable global state, it is best to set it once at the start of a session and not change it, or to setup and solve a problem within a single

```julia
setprecision(2000) do
    ...
end
```

block, to avoid any potentially confusing behavior caused by mixing precisions.
