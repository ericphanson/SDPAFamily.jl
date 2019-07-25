# SDPA_GMP

[![Build Status](https://travis-ci.com/ericphanson/SDPA_GMP.jl.svg?branch=master)](https://travis-ci.com/ericphanson/SDPA_GMP.jl)
[![Codecov](https://codecov.io/gh/ericphanson/SDPA_GMP.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/ericphanson/SDPA_GMP.jl)

A first attempt at using SDPA-GMP in Julia.

The binaries are built using BinaryProvider at the repo <https://github.com/ericphanson/SDPA_GMP_Builder>.

The plan is to follow [SDPA.jl](https://github.com/JuliaOpt/SDPA.jl) and use [CxxWrap.jl](https://github.com/JuliaInterop/CxxWrap.jl) to wrap SDPA-GMP for use in Julia, and write a [MathOptInterface.jl](https://github.com/JuliaOpt/MathOptInterface.jl) interface.

Currently, the library is not wrapped, so instead this package (`SDPA_GMP.jl`) just provides a binary installation and shorthand access to it via the Julia function `sdpa_gmp`, e.g.
```
sdpa_gmp_binary(`-o example1.result -dd example1.dat`)
```
following the example of [FFMPEG.jl](https://github.com/JuliaIO/FFMPEG.jl).
