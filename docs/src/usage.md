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
* `presolve`: whether or not to run a presolve routine to eliminate equality
  constraints and remove linearly dependent constraints. See below for more
  details. Defaults to `true`. Note,
  `presolve=true` is required to pass many of the tests for this package;
  linearly independent constraints are an assumption of the SDPA-family, but
  constraints generated from high level modelling languages often do have linear
  dependence between them.
* `presolve_tol`: how exactly a constraint must be a combination of others for
  the presolve to treat it as redundant, relative to the sizes of the
  coefficients. Defaults to `64 * eps(T)`; smaller values have no effect. The
  problem data are taken as exact numbers of type `T`, so if e.g. `Float64` data
  are used in a `BigFloat` model, `presolve_tol = 64 * eps(Float64)` may be
  needed for redundant constraints to be detected.
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

## Using presolve

SDPA-GMP will emit `cholesky miss condition :: not positive definite` error when
the problem data contain linearly dependent constraints. By default
(`presolve=true`), we use a presolve routine which detects such constraints by
Gaussian elimination, and which also uses the equality constraints to eliminate
variables. The redundant constraints are omitted from the problem formulation
and the corresponding redundant decision variables are set to 0 in the final
result. The time taken to perform the presolve step depends on the size of the
problem. The elimination is done in extended precision, and the solution
returned by the solver is checked against the original problem.

If that check fails, the termination status is `NUMERICAL_ERROR`; this usually
means the problem is badly scaled, since the tolerances of the SDPA solvers are
absolute. If the solver stops at its iteration limit with feasible points whose
objective values do not agree, the status is `ITERATION_LIMIT`. If the presolve
cannot confirm its own reduction, it emits a warning and passes the problem to
the solver without the reduction, as with `presolve=false`.

```@setup presolve
using SDPAFamily, Test, Convex, LinearAlgebra
```

We demonstrate `presolve` using a small problem in which a matrix variable `A`
only enters through its symmetric part `A + A'`. The optimal value is 1. When
`presolve` is disabled,
SDPA solvers will terminate prematurely due to linear dependence in the input
constraints. Note, however, that this does not necessarily happen. Empirically,
for our test cases, solvers' intolerance to redundant constraints increases from
`:sdpa` to `:sdpa_gmp`.

```@repl presolve
A = Variable(2, 2);
problem = minimize(tr(A), [A + A' ⪰ Diagonal(ones(2))]; numeric_type = BigFloat);
solve!(problem, () -> SDPAFamily.Optimizer(presolve = false))
```

Applying presolve helps by removing the redundant variable from the final input
file.

```@repl presolve
solve!(problem, () -> SDPAFamily.Optimizer(presolve = true))
@test problem.optval ≈ 1 atol=1e-30
```

We see we have recovered the true answer to a tolerance of $10^{-30}$.

## Changing parameters & solving at very high precision

```@setup convexquantum
using SDPAFamily, Test, Convex
```

Returning to the problem defined in [Optimal guessing probability for a pair of
quantum states](@ref), we can also increase the precision by changing the default parameters:

```@repl convexquantum
opt = () -> SDPAFamily.Optimizer(
    presolve = true,
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
