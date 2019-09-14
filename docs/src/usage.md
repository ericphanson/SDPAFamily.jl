# Usage

The main object of interest supplied by this package is `SDPAFamily.Optimizer{T}()`. Here, `T` is a numeric type which defaults to `BigFloat`. Keyword arguments may also be passed to `SDPAFamily.Optimizer()`:

* `variant`: either `:sdpa_gmp`, `:sdpa_qd`, or `:sdpa_dd`, to use SDPA-GMP, SDPA-QD, or SDPA-DD, respectively. Defaults to `:sdpa_gmp`.
* `silent`: a boolean to indicate whether or not to print output. Defaults to `true`.
* `presolve`: whether or not to run a presolve routine to remove linearly dependent constraints. See below for more details. Defaults to `false`. Note, `presolve=true` is required to pass many of the tests for this package; linearly independent constraints are an assumption of the SDPA-family, but constraints generated from high level modelling languages often do have linear dependence between them.
* `binary_path`: a string representing a path to the SDPA-GMP binary to use. The default is chosen at `build` time. To change the default binaries, see [Custom binary](@ref).
* `params_path`: a string representing a path to a parameter file named `param.sdpa` in the same folder as the binary if present. For details please refer to the [SDPA users manual](https://sourceforge.net/projects/sdpa/files/sdpa/sdpa.7.1.1.manual.20080618.pdf). The default parameters used by `SDPAFamily.jl` are here: <https://github.com/ericphanson/SDPAFamily.jl/blob/master/deps/>.

`SDPAFamily.Optimizer()` also accepts `variant = :sdpa` to use the non-high-precision SDPA binary. For general usage of the SDPA solver, use [SDPA.jl](https://github.com/JuliaOpt/SDPA.jl) which uses the C++ library to interface directly with the SDPA binary.

### Using a number type other than `BigFloat`

`SDPAFamily.Optimizer()` uses `BigFloat` for problem data and solution by default. To use, for example, `Float64` instead, simply call `SDPAFamily.Optimizer{Float64}()`. However, this may cause underflow errors when reading the solution file, particularly when `MathOptInterface` bridges are used. Note that with `MathOptInterface` all the problem data must be parametrised by the same number type, i.e. `Float64` in this case.

### Using presolve

SDPA-GMP will emit `cholesky miss condition :: not positive definite` error when the problem data contain linearly dependent constraints. Set `presolve=true` to use a presolver which tries to detect such constraints by Gaussian elimination. The redundant constraints are omitted from the problem formulation and the corresponding redundant decision variables are set to 0 in the final result. The time taken to perform the presolve step depends on the size of the problem. This process is numerically unstable, so is disabled by default. It does help quite a lot with some problems, however. 

```@setup convexquantum
using SDPAFamily, Test
using Pkg
Pkg.add(PackageSpec(name="Convex", url="https://github.com/ericphanson/Convex.jl", rev="MathOptInterface"));
using Convex
E12, E21 = ComplexVariable(2, 2), ComplexVariable(2, 2);
s1, s2 = [big"0.25" big"-0.25"*im; big"0.25"*im big"0.25"], [big"0.5" big"0.0"; big"0.0" big"0.0"];
p = Problem{BigFloat}(:minimize, real(tr(E12 * (s1 + 2 * s2) + E21 * (s2 + 2 * s1))), [E12 ⪰ 0, E21 ⪰ 0, E12 + E21 == Diagonal(ones(2)) ]);
opt = SDPAFamily.Optimizer(presolve = true, variant = :sdpa_gmp, silent = true);

```
We demonstrate `presolve` using the problem defined in [Examples](@ref). When `presolve` is disabled, SDPA solvers will terminate prematurely due to linear dependence in the input constraints. Note, however, that this does not necessarily happen. Empirically, for our test cases, solvers' intolerance increases from `:sdpa` to `:sdpa_gmp`.
```@repl convexquantum
opt = SDPAFamily.Optimizer(presolve = false)
solve!(p, SDPAFamily.Optimizer())
```
Applying presolve helps by removing 8 redundant constraints from the final input file.
```@repl convexquantum
opt.presolve = true; opt.silent = true;
solve!(p, opt)
@test p.optval ≈ (6 - sqrt(big"2.0"))/4 atol=1e-30
SDPAFamily.presolve(opt)
```