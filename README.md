# SDPA_GMP

[![Build Status](https://travis-ci.com/ericphanson/SDPA_GMP.jl.svg?branch=master)](https://travis-ci.com/ericphanson/SDPA_GMP.jl)
[![Codecov](https://codecov.io/gh/ericphanson/SDPA_GMP.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/ericphanson/SDPA_GMP.jl)

An interface to using SDPA-GMP, SDPA-DD, and SDPA-QD in Julia (<http://sdpa.sourceforge.net>). Call `SDPA_GMP.Optimizer()` to use this wrapper via `MathOptInterface`, which is an intermediate layer between low-level solvers (such as SDPA-GMP, SDPA-QD, and SDPA-DD) and high level modelling languages, such as [JuMP.jl](https://github.com/JuliaOpt/JuMP.jl) and [Convex.jl](https://github.com/JuliaOpt/Convex.jl/).

JuMP currently only supports `Float64` numeric types, which means that problems can only be specified to 64-bits of precision, and results can only be recovered at that level of precision, when using JuMP. This is tracked in the issue [JuMP#2025](https://github.com/JuliaOpt/JuMP.jl/issues/2025).

Convex.jl does not yet officially support MathOptInterface; this issue is tracked at [Convex.jl#262](https://github.com/JuliaOpt/Convex.jl/issues/262). However, there is a work-in-progress branch which can be added to your Julia environment via
```julia
] add https://github.com/ericphanson/Convex.jl#MathOptInterface
```
which can be used to solve problems with the solvers from this package.

## Installation

This package is not yet registered in the General registry. To install, type `]` in the julia command prompt, then execute

```julia
pkg> add https://github.com/ericphanson/SDPA_GMP.jl
```

### Automatic binary installation

If you are on MacOS or Linux, this package will attempt to automatically download the SDPA-GMP, SDPA-DD, and SDPA-QD binaries, built by [SDPA_GMP_Builder.jl](https://github.com/ericphanson/SDPA_GMP_Builder). The SDPA-GMP binary is patched from the official SDPA-GMP source to allow printing more digits, in order to recover high-precision output.

SDPA-GMP does not compile on Windows. However, it can be used via the [Windows Subsystem for Linux](https://docs.microsoft.com/en-us/windows/wsl/about)(WSL). If you have WSL installed, then SDPA_GMP.jl will try to automatically detect this and use an appropriate binary, called via WSL. This binary can be found at the repo <https://github.com/ericphanson/SDPA_GMP_Builder>, and is built on WSL from the source code at <https://github.com/ericphanson/sdpa-gmp>.

SDPA-{GMP, QD, DD} are each available under a GPLv2 license, which can be found here: <https://github.com/ericphanson/SDPA_GMP_Builder/blob/master/deps/COPYING>.

### Custom binary

If you would like to use a different binary, set the enviromental variable `JULIA_SDPA_GMP_PATH` to the folder containing the binary you would like to use, and then build the package. This can be done in Julia by, e.g.,

```julia
ENV["JULIA_SDPA_GMP_PATH"] = "/path/to/folder/"
import Pkg
Pkg.build("SDPA_GMP")
```

and that will configure this package to use that binary by default. If your custom location is via WSL on Windows, then also set `ENV["JULIA_SDPA_GMP_WSL"] = "TRUE"` so that SDPA_GMP.jl knows to adjust the paths to the right format.

It is recommended to patch SDPA-GMP (as was done in <https://github.com/ericphanson/sdpa-gmp>) in order to allow printing more digits.

* For source code downloaded from the official website (dated 20150320), modify the `P_FORMAT` string at line 23 in `sdpa_struct.h` so that the output has a precision no less than 200 bits (default) or precision specified by the parameter file. 
* For source code downloaded from its [GitHub repository](https://github.com/nakatamaho/sdpa-gmp), specify the print format string in `param.sdpa` as described in the [SDPA users manual](https://sourceforge.net/projects/sdpa/files/sdpa/sdpa.7.1.1.manual.20080618.pdf).

Other information about compiling SDPA-GMP binary can be found [here](https://sourceforge.net/projects/sdpa/files/sdpa-gmp/sdpa-gmp.7.1.2-install.txt).

## Usage

The main object of interest supplied by this package is `SDPA_GMP.Optimizer{T}()`. Here, `T` is a numeric type which defaults to `BigFloat`. Keyword arguments may also be passed to `SDPA_GMP.Optimizer`:

* `variant`: either `:gmp`, `:qd`, or `:dd`, to use SDPA-GMP, SDPA-QD, or SDPA-DD, respectively. Defaults to `:gmp`.
* `silent`: a boolean to indicate whether or not to print output. Defaults to `true`.
* `presolve`: whether or not to run a presolve routine to remove linearly dependent constraints. See below for more details. Defaults to `false`.
* `binary_path`: a string representing a path to the SDPA-GMP binary to use. The default is chosen at `build` time.
* `params_path`: a representing a path to a parameter file named `param.sdpa` in the same folder as the binary if present. For details please refer to the [SDPA users manual](https://sourceforge.net/projects/sdpa/files/sdpa/sdpa.7.1.1.manual.20080618.pdf). The default parameters used by SDPA_GMP.jl are here: <https://github.com/ericphanson/SDPA_GMP.jl/blob/master/deps/param.sdpa>.

### Example

[todo]

### Using a number type other than `BigFloat`

SDPA-GMP.jl uses `BigFloat` for problem data and solution by default. To use, for example, `Float64` instead, simply call `SDPA_GMP.Optimizer{Float64}()`. However, this may cause underflow errors when reading the solution file, particularly when `MathOptInterface` bridges are used. Note that with `MathOptInterface` all the problem data must be parametrised by the same number type, i.e. `Float64` in this case.

### Using presolve

SDPA-GMP will emit `cholesky miss condition :: not positive definite` error if the problem data contain linearly dependent constraints. Set `presolve=true`, to use a presolver to try to detect such constraints by Gaussian elimination. The redundant constraints are omitted from the problem formulation and the corresponding decision variables are set to 0 in the final result. The time taken to perform the presolve step depends on the size of the problem. This process is numerically unstable, so is disabled by default. It does help quite a lot with some problems, however.
