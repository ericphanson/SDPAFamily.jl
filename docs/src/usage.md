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
* `presolve`: whether or not to run a presolve routine to remove linearly
  dependent constraints. See below for more details. Defaults to `false`. Note,
  `presolve=true` is required to pass many of the tests for this package;
  linearly independent constraints are an assumption of the SDPA-family, but
  constraints generated from high level modelling languages often do have linear
  dependence between them.
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
the problem data contain linearly dependent constraints. Set `presolve=true` to
use a presolve routine which tries to detect such constraints by Gaussian
elimination. The redundant constraints are omitted from the problem formulation
and the corresponding redundant decision variables are set to 0 in the final
result. The time taken to perform the presolve step depends on the size of the
problem. This process is numerically unstable, so is disabled by default. It
does help quite a lot with some problems, however.

```@setup convexquantum
using SDPAFamily, Test
using Pkg
Pkg.add(PackageSpec(name="Convex", url="https://github.com/ericphanson/Convex.jl", rev="MathOptInterface"));
using Convex
ρ₁ = Complex{BigFloat}[1 0; 0 0]
ρ₂ = (1//2)*Complex{BigFloat}[1 -im; im 1]
E₁ = ComplexVariable(2, 2);
E₂ = ComplexVariable(2, 2);
problem = maximize( real((1//2)*tr(ρ₁*E₁) + (1//2)*tr(ρ₂*E₂)),
                    [E₁ ⪰ 0, E₂ ⪰ 0, E₁ + E₂ == Diagonal(ones(2))];
                    numeric_type = BigFloat );
opt = SDPAFamily.Optimizer(presolve = true, variant = :sdpa_gmp, silent = true);
```

We demonstrate `presolve` using the problem defined in [Optimal guessing
probability for a pair of quantum states](@ref). When `presolve` is disabled,
SDPA solvers will terminate prematurely due to linear dependence in the input
constraints. Note, however, that this does not necessarily happen. Empirically,
for our test cases, solvers' intolerance to redundant constraints increases from
`:sdpa` to `:sdpa_gmp`.

```@repl convexquantum
solve!(problem, SDPAFamily.Optimizer(presolve = false))
```

Applying presolve helps by removing 8 redundant constraints from the final input
file.

```@repl convexquantum
solve!(problem, SDPAFamily.Optimizer(presolve = true))
@test problem.optval ≈ 1//2 + 1/(2*sqrt(big(2))) atol=1e-30
```

We see we have recovered the true answer to a tolerance of $10^{-30}$.

## Changing parameters & solving at very high precision

Continuing the above example, we can also increase the precision by changing the default parameters:

```@repl convexquantum
opt = SDPAFamily.Optimizer(
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

In this case, since `1` and `1/2` are exactly representable by floating point numbers, it is enough to specify `setprecision(2000)` before the `solve!` call (so the `ρ₁` and `ρ₂` are only constructed at 256 bits of high precision), but it is good practice to set the precision at the start for the whole problem. Moreover, since the precision is mutable global state, it is best to set it once at the start of a session and not change it, to avoid any potentially confusing behavior, or setup and solve problems within a single

```julia
setprecision(2000) do
    ...
end
```

block.
