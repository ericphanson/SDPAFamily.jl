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

SDPA solves the pair of problems

```math
\begin{aligned}
(P) \quad & \text{minimize } c^T x && \text{subject to } F_1 x_1 + \dots + F_m x_m - F_0 \succeq 0, \\
(D) \quad & \text{maximize } F_0 \bullet Y && \text{subject to } F_i \bullet Y = c_i \text{ for } i = 1, \dots, m, \ Y \succeq 0.
\end{aligned}
```

The variables of (P) are the variables of the MathOptInterface model, and each
constraint of the model is a block of the matrices ``F_i``. Models built with
Convex.jl, JuMP or MathOptInterface bridges often have two properties that the
primal-dual interior point method of SDPA handles poorly:

1. Equality constraints. Without a presolve, MathOptInterface writes each
   equality constraint ``a^T x + b = 0`` as the two inequalities
   ``a^T x + b \ge 0`` and ``-a^T x - b \ge 0``. Then (P) has no strictly
   feasible point.
2. Linearly dependent matrices ``F_i``. Then the equality constraints of (D) are
   linearly dependent, and the Schur complement matrix that SDPA factorizes in
   each iteration is singular. SDPA then prints `cholesky miss condition`
   warnings, stops at the iteration limit, or returns an inaccurate point.

With `presolve = true` (the default), SDPAFamily removes both before it writes
the input file. The computations use the number type `T` of
`SDPAFamily.Optimizer{T}`, so with `BigFloat` they are carried out at the
precision of `BigFloat`.

1. Gaussian elimination on the equality constraints expresses one variable for
   each linearly independent equality constraint as an affine function of the
   other variables. These expressions are substituted into the other
   constraints and the objective. A linearly dependent equality constraint is
   dropped if it is consistent with the others. If it is not, the model is
   infeasible: the termination status is `INFEASIBLE`, the dual status is
   `INFEASIBILITY_CERTIFICATE`, and SDPA is not called.
2. Gaussian elimination on the matrices ``F_i`` of the remaining problem finds a
   maximal linearly independent subset. The other variables are fixed to zero.
   If ``F_k = \sum_i a_i F_i`` but ``c_k \ne \sum_i a_i c_i``, then (D) is
   infeasible: the termination status is `DUAL_INFEASIBLE`, the primal status
   is `INFEASIBILITY_CERTIFICATE` (an unbounded ray), and SDPA is not called.

After SDPA has solved the reduced problem, the solution is mapped back to the
model. The values of the eliminated variables are computed from the other
variables. The primal and dual values of the cone constraints need no
correction. The dual values ``y_E`` of the equality constraints are computed from
the dual values ``y_K`` of the cone constraints by solving
``A_E^T y_E = c - A_K^T y_K``, where ``A_E`` and ``A_K`` are the coefficient
matrices of the equality and cone constraints. Dropped equality constraints
have the dual value zero.

If the equality constraints determine all variables, SDPA is not called. The
termination status is `OPTIMAL` if the cone constraints hold up to a relative
tolerance of `sqrt(eps(T))`, `ALMOST_OPTIMAL` (with primal status
`NEARLY_FEASIBLE_POINT`) if they hold up to `sqrt(eps(Float64))`, and
`INFEASIBLE` (without a certificate) otherwise.

A row or column counts as linearly dependent if, after elimination, its largest
remaining entry is at most `sqrt(eps(T))` times its largest original entry. Rows
or columns that are close to dependent but above this tolerance are kept.

With `presolve = false`, the problem is passed to SDPA unchanged, with equality
constraints bridged into pairs of inequalities.

The following problem has two variables whose matrices ``F_i`` are equal:

```@repl presolve
using SDPAFamily, Convex
x = Variable(); y = Variable();
problem = minimize(x + y, [x + y >= 1]; numeric_type = BigFloat);
solve!(problem, () -> SDPAFamily.Optimizer(presolve = false, silent = true))
problem.status
solve!(problem, () -> SDPAFamily.Optimizer(presolve = true, silent = true))
problem.status
problem.optval
```

Without the presolve, SDPA-GMP stops at the iteration limit. With the presolve,
one of the two variables is fixed to zero and SDPA-GMP finds the optimal value
1.

To see what the presolve removed, use `verbose = SDPAFamily.VERBOSE`, which
prints the number of eliminated and fixed variables.

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
