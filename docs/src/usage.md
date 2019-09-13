# Usage

The main object of interest supplied by this package is `SDPAFamily.Optimizer{T}()`. Here, `T` is a numeric type which defaults to `BigFloat`. Keyword arguments may also be passed to `SDPAFamily.Optimizer`:

* `variant`: either `:gmp`, `:qd`, or `:dd`, to use SDPA-GMP, SDPA-QD, or SDPA-DD, respectively. Defaults to `:gmp`.
* `silent`: a boolean to indicate whether or not to print output. Defaults to `true`.
* `presolve`: whether or not to run a presolve routine to remove linearly dependent constraints. See below for more details. Defaults to `false`. Note, `presolve=true` is required to pass many of the tests for this package; linearly independent constraints are an assumption of the SDPA-family, but constraints generated from high level modelling languages often do have linear dependence between them.
* `binary_path`: a string representing a path to the SDPA-GMP binary to use. The default is chosen at `build` time.
* `params_path`: a string representing a path to a parameter file named `param.sdpa` in the same folder as the binary if present. For details please refer to the [SDPA users manual](https://sourceforge.net/projects/sdpa/files/sdpa/sdpa.7.1.1.manual.20080618.pdf). The default parameters used by SDPAFamily.jl are here: <https://github.com/ericphanson/SDPAFamily.jl/blob/master/deps/param.sdpa>.

### Using a number type other than `BigFloat`

SDPA-GMP.jl uses `BigFloat` for problem data and solution by default. To use, for example, `Float64` instead, simply call `SDPAFamily.Optimizer{Float64}()`. However, this may cause underflow errors when reading the solution file, particularly when `MathOptInterface` bridges are used. Note that with `MathOptInterface` all the problem data must be parametrised by the same number type, i.e. `Float64` in this case.

### Using presolve

SDPA-GMP will emit `cholesky miss condition :: not positive definite` error if the problem data contain linearly dependent constraints. Set `presolve=true`, to use a presolver to try to detect such constraints by Gaussian elimination. The redundant constraints are omitted from the problem formulation and the corresponding decision variables are set to 0 in the final result. The time taken to perform the presolve step depends on the size of the problem. This process is numerically unstable, so is disabled by default. It does help quite a lot with some problems, however.
