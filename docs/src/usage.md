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
the input file:

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

### Arithmetic and tolerances

The problem data of type `T` are taken as exact numbers and converted exactly to
a working type with more precision: `Float64x3` from
[MultiFloats.jl](https://github.com/dzhang314/MultiFloats.jl) (about 157 bits)
for `T = Float64`, and `BigFloat` with `2 * precision(T) + 64` bits otherwise.
The presolve computes in the working type. In this section, ``\epsilon_T`` is the
machine epsilon of `T` (`eps(T)`, for `BigFloat` at the precision set by the
caller), ``\epsilon_W`` that of the working type, and
``\tau = \sqrt{\epsilon_T \epsilon_W}``.

Each decision compares a computed combination ``\sum_j w_j a_j`` of data vectors
``a_j`` (rows of the equality constraints, or columns of the constraint matrix)
with the final coefficients ``w_j`` against the same data. An entry ``r`` of the
combination counts as zero at tolerance ``t`` if
``|\sum_j w_j a_{jr}| \le t \sum_j |w_j| |a_{jr}| + \tau \sum_j |w_j| \max_s |a_{js}|``.
The second term covers rounding errors of the working type that spread across
the entries of a combination. There are two tolerances:

* ``\tau`` (rounding): entries that count as zero at tolerance ``\tau`` consist
  of rounding errors of the working type. They are removed during the
  elimination and the substitution.
* ``\delta = 64 \epsilon_D`` (data): a row of the equality constraints, or a
  matrix ``F_k``, is linearly dependent on the others if all its entries count as
  zero at tolerance ``\delta`` after the elimination (confirmed by recomputing
  the combination from the data), and a dependent row or matrix is consistent if
  the same combination of the right-hand sides (objective coefficients) counts as
  zero at tolerance ``\delta``. Here ``\epsilon_D = \epsilon_T``, except that
  ``\epsilon_D = \epsilon_{\text{Float64}}`` if `T` has more precision than
  `Float64` and all coefficients of the model are `Float64` numbers: such data
  were most likely computed in `Float64` (for example decimal coefficients like
  `0.1` widened to `BigFloat`), and relations between them hold only up to
  rounding to `Float64`.

### Scaling

Before the eliminations, each row of the equality and cone constraints is scaled
to a largest coefficient near 1, so that the sizes of combinations of rows are
comparable. After the reduction, the problem passed to SDPA is equilibrated by
Ruiz's algorithm, which divides the rows (with the right-hand side counted in the
row norms) and the variables repeatedly by the square roots of their largest
coefficients; the rows of a positive semidefinite block share one factor, so
that the scaling keeps the cone. Then the right-hand side and the objective are
scaled so that their largest entries are between 1 and 2. All factors are powers
of two, so the scaled data are exact. The solution is scaled back before the
checks below.

### Checks and fallback

Every result that the presolve returns without calling SDPA is checked against
the data before it is returned, at the tolerance ``\delta`` of the decisions
that led to it: a certificate of infeasibility must satisfy ``A^T y = 0`` and
``b^T y < 0``, a ray must satisfy ``A_E d = 0``, ``A_K d \in K`` and
``c^T d < 0``, and the expressions for the eliminated variables are checked
against the original equality constraints.

If no variable is left (because the equality constraints determine all
variables, or because the model has no variables with nonzero coefficients),
SDPA is not called, and the presolve checks the constraints at the only
remaining point, computed from the data. A `Nonnegatives` entry with value
``b_i`` and magnitude ``m_i`` (the sum of the absolute values of the terms that
give ``b_i``) is accepted if ``b_i \ge -\epsilon m_i``, and a positive
semidefinite block ``B`` if ``B + \epsilon \|M\|_F I`` is positive
semidefinite, where ``M`` holds the magnitudes of the entries of ``B``. Here
``\epsilon`` is the larger of ``\delta`` and the constraint tolerance
`epsilonStar` of the solver parameters (zero when the parameters are given as a
file or as `UNSTABLE_BUT_FAST` or `STABLE_BUT_SLOW`). Semidefiniteness is tested
by symmetric elimination that accepts zero pivots with zero rows, so the zero
matrix and other singular positive semidefinite matrices are accepted. If all
constraints are accepted, the termination status is `OPTIMAL`. Otherwise it is
`INFEASIBLE`; if a `Nonnegatives` entry is violated, the dual status is
`INFEASIBILITY_CERTIFICATE`, and otherwise `NO_SOLUTION`.

After SDPA has solved the reduced problem, the solution is mapped back and
checked:

* The equality constraints that were kept must hold up to ``64 \epsilon_T``
  relative to the sizes of their terms (the rounding of the solution to `T`).
  A constraint that was dropped as dependent must hold up to ``2\delta`` times
  the terms of the combination it was found to depend on, plus the residuals of
  the constraints in that combination.
* When the dual status is `FEASIBLE_POINT`, the relative residual of the dual
  equations ``A^T y = c`` of the model must not exceed 100 times the relative
  residual of the dual equations of the reduced problem plus ``\delta``. The
  mapping back is exact up to rounding, so this only fails when the elimination
  amplifies the solver's error, for example for nearly dependent constraints.

If any check fails, the presolve falls back: the problem is passed to SDPA
without reduction and without scaling, with each equality constraint written as
two inequalities, and a warning is printed unless `silent = true`. SDPA may then
fail on the problem, but the result is the solver's own.

### Cost

Eliminating variables can fill in sparse rows of the constraints. The presolve
chooses as pivots variables with few nonzero coefficients, among those whose
coefficient is at least a tenth of the largest one in the row. In a test with
400 variables, 200 random equality constraints with 3 nonzero coefficients
each, and `x ≥ 0`, the SDPA input file was 84 KiB with the presolve and 34 KiB
without it. For 1000 random sparse equality constraints on a 60×60 matrix
variable, the presolve took about 0.9 seconds with `Float64` and 16 seconds with
`BigFloat`.

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
