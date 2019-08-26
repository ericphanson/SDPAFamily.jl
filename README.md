# SDPA_GMP

[![Build Status](https://travis-ci.com/ericphanson/SDPA_GMP.jl.svg?branch=master)](https://travis-ci.com/ericphanson/SDPA_GMP.jl)
[![Codecov](https://codecov.io/gh/ericphanson/SDPA_GMP.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/ericphanson/SDPA_GMP.jl)

A first attempt at using [SDPA-GMP](http://sdpa.sourceforge.net/download.html#sdpa-gmp) in Julia. Call `SDPA_GMP.Optimizer()` to use this wrapper via `MathOptInterface`. 

## Installation

This package is not yet registered in `METADATA.jl`. To install, type `]` in the julia command prompt, then execute

```
(v1.1) pkg> add https://github.com/ericphanson/SDPA_GMP.jl.git
```

After installing the package, make sure that the `sdpa_gmp` binary is in system PATH. Alternatively, modify the `sdpa_gmp_path` constant in `SDPA_GMP.jl` to a path that points to the `sdpa_gmp` binary. Further information about compiling SDPA-GMP binary can be found [here](https://sourceforge.net/projects/sdpa/files/sdpa-gmp/sdpa-gmp.7.1.2-install.txt). 

Note that before compilation it is often necessary to modify the SDPA-GMP source code so that the output values have sufficient precision. To do this, modify the `P_FORMAT` string at line 23 in `sdpa_struct.h` so that the output has a precision no less than the default precision (200 bits) or precision specified by the parameter file. 

## Usage

### Parameters

SDPA-GMP uses the parameter file named `param.sdpa` in the same folder as the binary if present. For details please refer to the [SDPA users manual](https://sourceforge.net/projects/sdpa/files/sdpa/sdpa.7.1.1.manual.20080618.pdf).

### Using a number type other than `Float64`

SDPA-GMP.jl uses `Float64` for problem data and solution by default. To use, for example, `BigFloat` instead, simply call `SDPA_GMP.Optimizer{BigFloat}()`. Note that with `MathOptInterface` all the problem data must be of the same type, i.e. `BigFloat` in this case. 

### Using presolve

SDPA-GMP will emit `cholesky miss condition :: not positive definite` error if the problem data contain linearly dependent constraints. By default, a presolver will try to detect such constraints by Gaussian elimination. The redundant constraints are omitted from the problem formulation and the corresponding decision variables are set to 0 in the final result. The time taken to perform the presolve step depends on the size of the problem. To disable, call `SDPA_GMP.Optimizer(presolve = false)` instead. 
