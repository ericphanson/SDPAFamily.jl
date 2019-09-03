# SDPA_GMP

[![Build Status](https://travis-ci.com/ericphanson/SDPA_GMP.jl.svg?branch=master)](https://travis-ci.com/ericphanson/SDPA_GMP.jl)
[![Codecov](https://codecov.io/gh/ericphanson/SDPA_GMP.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/ericphanson/SDPA_GMP.jl)

A first attempt at using [SDPA-GMP](http://sdpa.sourceforge.net/download.html#sdpa-gmp) in Julia. Call `SDPA_GMP.Optimizer()` to use this wrapper via `MathOptInterface`. 

## Installation

This package is not yet registered in the General registry. To install, type `]` in the julia command prompt, then execute

```julia
pkg> add https://github.com/ericphanson/SDPA_GMP.jl.git
```

### MacOS or Linux

If you are on MacOS or Linux, this package will attempt to automatically download the SDPA-GMP binary, built from <https://github.com/ericphanson/sdpa-gmp> (built by [SDPA_GMP_Builder.jl](https://github.com/ericphanson/SDPA_GMP_Builder)). This is patched from the official SDPA-GMP source to allow printing more digits, in order to recover high-precision output.

### Windows

SDPA-GMP does not compile on Windows. However, it can be used via the [Windows Subsystem for Linux](https://docs.microsoft.com/en-us/windows/wsl/about)(WSL). You can install the SDPA-GMP binary on your WSL installation by issuing the following commands in the WSL prompt:

```bash
mkdir SDPA_GMP
cd SDPA_GMP
wget https://github.com/ericphanson/SDPA_GMP_Builder/releases/download/v7.1.3/SDPA_GMP_Builder.v7.1.3.x86_64-linux-gnu.tar.gz
tar -xf SDPA_GMP_Builder.v7.1.3.x86_64-linux-gnu.tar.gz
cp usr/bin/sdpa_gmp /usr/bin/local/sdpa_gmp
```

If you have WSL installed and SDPA-GMP installed on your WSL path (as the above commands try to do), then SDPA_GMP.jl will try to automatically detect this and use that SDPA-GMP binary.

### Custom binary location

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

* `silent`: a boolean to indicate whether or not to print output. Defaults to `true`.
* `presolve`: whether or not to run a presolve routine to remove linearly dependent constraints. See below for more details. Defaults to `true`.
* `binary_path`: a string representing a path to the SDPA-GMP binary to use. The default is chosen at `build` time.
* `params_path`: a representing a path to a parameter file named `param.sdpa` in the same folder as the binary if present. For details please refer to the [SDPA users manual](https://sourceforge.net/projects/sdpa/files/sdpa/sdpa.7.1.1.manual.20080618.pdf). The default parameters used by SDPA_GMP.jl are here: <https://github.com/ericphanson/SDPA_GMP.jl/blob/master/deps/param.sdpa>.

### Example

[todo]

### Using a number type other than `BigFloat`

SDPA-GMP.jl uses `BigFloat` for problem data and solution by default. To use, for example, `Float64` instead, simply call `SDPA_GMP.Optimizer{Float64}()`. However, this may cause underflow errors when reading the solution file, particularly when `MathOptInterface` bridges are used. Note that with `MathOptInterface` all the problem data must be parametrised by the same number type, i.e. `Float64` in this case.

### Using presolve

SDPA-GMP will emit `cholesky miss condition :: not positive definite` error if the problem data contain linearly dependent constraints. By default, a presolver will try to detect such constraints by Gaussian elimination. The redundant constraints are omitted from the problem formulation and the corresponding decision variables are set to 0 in the final result. The time taken to perform the presolve step depends on the size of the problem. To disable, call `SDPA_GMP.Optimizer(presolve = false)` instead.
