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
* `presolve_tol`: the tolerance of the presolve for linear dependence and
  consistency, see [Arithmetic and tolerances](@ref). Defaults to `nothing`,
  which means `64 * eps(T)`.
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
`BigFloat` with twice the precision of `T` plus 64 bits (170 bits for
`T = Float64`); this is the working precision, in which the presolve computes.
In this section, ``p`` is the precision of `T` in bits and ``\epsilon_T`` its
machine epsilon (`eps(T)`, for `BigFloat` at the precision set by the caller).

There are two tolerances:

* Rounding residue. An exact cancellation of two numbers leaves a residue near
  the working precision. A computed number is set to zero if its terms cancelled
  so far that its leading bit lies more than ``1.5 p + 32`` bits below the
  leading bit of its largest term (111 bits for `Float64`). Numbers that follow
  exactly from the data are kept: a coefficient ``10^{-100}`` next to a
  coefficient 1 in the same row, or the coefficient ``s`` of ``L x - L y + s w``
  after the substitution ``x = y = w`` with ``s = 10^{-8}`` and ``L = 10^8``.
* Dependence, with ``\delta = 64 \epsilon_T`` by default. A row of the equality
  constraints (or a matrix ``F_k``) is linearly dependent on the rows (matrices)
  before it if a combination with coefficients ``w_j`` of the data vectors
  ``a_j`` satisfies ``|\sum_j w_j a_j| \le \delta \sum_j |w_j| |a_j|`` in every
  entry. The combination is found by the elimination and the condition is then
  evaluated from the data. It means that relative changes of the data of at most
  ``\delta`` make the vectors exactly dependent (the componentwise backward
  error of Oettli and Prager), and it does not change when rows or variables
  are scaled. A dependent row or matrix is consistent if the same combination of
  the right-hand sides (objective coefficients) is zero in the same sense,
  relative to the terms of the right-hand sides (objective coefficients) of the
  model.

The eliminations are sparse LU factorizations with threshold pivoting: a pivot is
at least a tenth of the largest entry of its row in the equilibrated model, and
among those the variable with the fewest nonzero coefficients is taken.

The data are taken as exact numbers of type `T`, also when they happen to be
`Float64` numbers: the model with constraints ``(10^8+1) x + 10^8 y = 1`` and
``10^8 x + (10^8-1) y = 0`` has the unique solution ``x = 1 - 10^8``, ``y = 10^8``
also in `BigFloat`. Data that were rounded to a coarser type, such as decimal
coefficients like `0.1` computed in `Float64` and then converted to `BigFloat`,
satisfy linear relations only up to that rounding. For such data, pass
`presolve_tol = 64 * eps(Float64)`; then rows and matrices that are dependent up
to the rounding are treated as dependent. A larger `presolve_tol` treats more
nearly dependent constraints as dependent, which changes the problem by up to
that relative amount.

### Scaling

The reduced problem is scaled before it is passed to SDPA, as conic solvers such
as SCS do. All factors are powers of two, so the scaled data are exact.

1. The constraint matrix is equilibrated by Ruiz's algorithm, which divides the
   rows and the variables repeatedly by the square roots of their largest
   coefficients. The rows of a positive semidefinite block share one factor, so
   that the scaling keeps the cone.
2. The right-hand side is normalized to a largest entry between 1 and 2, by
   multiplying the rows by a factor and dividing the variables by it. This is
   done separately for each group of rows and variables that are connected by
   nonzero coefficients.
3. The objective is multiplied by one factor ``\sigma`` so that its largest
   coefficient is between 1 and 2.

SDPA's tolerances are absolute, so they apply to the scaled problem. Its
feasibility tolerance is relative to the equilibrated rows and variables: a
variable is accurate relative to the unit that the scaling gave it, not relative
to its own value. Its gap tolerance `epsilonStar` is relative to
``\max(1, |\text{objective}|)`` of the scaled problem, that is, to
``\max(1/\sigma, |\text{objective}|)`` in the units of the model.

Rows of the cone constraints without coefficients (constants, for example from
a variable that Convex.jl fixed, or rows whose variables were eliminated) are
checked at the point given by the equality constraints, with the tolerance
``\epsilon`` described below. If they hold, their constants are replaced by 1
(or a constant positive semidefinite block by the identity matrix), which keeps
the feasible set and gives SDPA a strictly feasible row, and they do not take
part in the scaling. Their primal values are the original constants and their
dual values are zero. If a constant `Nonnegatives` row is violated, the
termination status is `INFEASIBLE` with a certificate, and SDPA is not called.

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
``\epsilon`` is the larger of ``\delta`` and the tolerance `epsilonStar` of the
solver parameters. Semidefiniteness is tested by symmetric elimination that
accepts zero pivots with zero rows, so the zero matrix and other singular
positive semidefinite matrices are accepted. If all constraints are accepted,
the termination status is `OPTIMAL`. Otherwise it is `INFEASIBLE`; if a
`Nonnegatives` entry is violated, the dual status is
`INFEASIBILITY_CERTIFICATE`, and otherwise `NO_SOLUTION`.

After SDPA has solved the reduced problem, its result is checked in two steps.

First, on the scaled problem that SDPA was given, with the tolerance
``\gamma = \sqrt{\max(\delta, \epsilon^*, \epsilon')}``, where ``\epsilon^*`` and
``\epsilon'`` are the tolerances `epsilonStar` and `epsilonDash` of the solver
parameters (``\gamma = 10^{-4}`` if they are given in a file that cannot be
read). This asks for half of the digits that SDPA was asked for; it catches
results that are wrong, not results that are inaccurate.

* A primal point, a dual point, and a result reported as optimal (SDPA's
  statuses `pdOPT` and `pdFEAS`) must satisfy the stopping criteria of SDPA with
  ``\gamma`` in place of its tolerances: the cone constraints, the dual
  equations ``A^T y = c``, and the duality gap relative to
  ``\max(1, |\text{objective}|)``.
* The gap that SDPA's tolerance admits must also be small in the units of the
  model: ``\epsilon^* \max(1, |\text{objective}|) / \sigma`` must be at most
  ``\gamma \max(1, |\text{objective of the model}|)``. This fails when the
  objective value is far smaller than the scaling could tell from the data,
  for example when a variable with a large cost can take large values on the
  feasible set but is zero at the optimum.
* SDPA does not return rays. It reports unboundedness when the objective of a
  feasible iterate passes its parameters `lowerBound` or `upperBound`
  (``\pm 10^5`` by default), and returns that iterate. A certificate of
  infeasibility ``y`` is accepted if ``b^T y < 0`` and
  ``\|A^T y\|_\infty \le 10^{-3} |b^T y|``, and a primal ray ``d`` if
  ``c^T d < 0`` and ``A d`` is in the cone up to ``10^{-3} |c^T d|``. This is the
  criterion of solvers that detect infeasibility from approximate certificates
  (SCS, OSQP, Clarabel). An accepted certificate excludes all solutions of the
  scaled problem with 1-norm below 1000. In the scaled problem the right-hand
  side and the objective are at most 2, so an iterate at SDPA's default bound
  has a ratio below ``2 \cdot 10^{-5}``.

Second, after the solution is mapped back, what the presolve removed is checked
at the returned solution, each relative to the sizes of its own terms:

* The equality constraints that were kept must hold up to ``64 \epsilon_T`` (the
  rounding of the solution to `T`). A constraint that was dropped as dependent
  must hold up to ``2\delta`` times the terms of the combination it was found
  to depend on, plus the residuals of the constraints in that combination.
* The dual equations of the eliminated variables must hold up to
  ``64 \epsilon_T``, and those of the variables that were fixed to zero up to
  ``2\delta`` times the terms of their combination, plus the residuals of the
  variables in that combination.

These checks fail when the elimination amplifies the solver's error, for
example for nearly dependent constraints.

If SDPA stopped at its iteration limit with feasible iterates (`pdFEAS`) and only
the gap criteria fail, the termination status is `ITERATION_LIMIT` and the primal
and dual statuses are `FEASIBLE_POINT`. If the first step fails otherwise, the
termination status is `NUMERICAL_ERROR`, the primal
and dual statuses are `UNKNOWN_RESULT_STATUS`, and `RawStatusString` gives the
status that SDPA reported and the check that failed. The same holds for a
solution that passes the first step but is not finite in `T` after it is scaled
back (for `T = Float64`, a variable beyond about ``1.8 \cdot 10^{308}``).

If the reduction itself is in doubt, the presolve falls back: the problem is
passed to SDPA without the reductions, with each equality constraint written as
two inequalities, scaled in the same way, and a warning is printed unless
`silent = true`. This happens if the working precision was not sufficient for
the eliminations, if a result of the presolve fails its validation, or if the
second step fails. The result of the fallback is checked on its scaled problem
as in the first step.

### Cost

Eliminating variables can fill in sparse rows of the constraints. The presolve
chooses as pivots variables with few nonzero coefficients, among those whose
coefficient is at least a tenth of the largest one in the row of the
equilibrated model. In a test with 400 variables, 200 random equality
constraints with 3 nonzero coefficients each, and `x ≥ 0`, the SDPA input file
was 86 KiB with the presolve and 34 KiB without it.

The time of the presolve is that of sparse Gaussian elimination in `BigFloat`
arithmetic. For 1000 random equality constraints with 5 nonzero coefficients on
a 60×60 matrix variable, and for 2000 random equality constraints with 3 nonzero
coefficients on 4000 nonnegative variables, it took about 0.1 seconds. For 400
dense equality constraints on 420 variables it took about 5 seconds with
`T = Float64` and 7 seconds with `T = BigFloat` at 256 bits.

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
