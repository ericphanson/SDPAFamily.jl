# SDPAFamily

[![Build Status](https://github.com/ericphanson/SDPAFamily.jl/workflows/CI/badge.svg)](https://github.com/ericphanson/SDPAFamily.jl/actions)
[![Codecov](https://codecov.io/gh/ericphanson/SDPAFamily.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/ericphanson/SDPAFamily.jl)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://ericphanson.github.io/SDPAFamily.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://ericphanson.github.io/SDPAFamily.jl/dev)

An interface to using SDPA-GMP, SDPA-DD, and SDPA-QD in Julia
(<http://sdpa.sourceforge.net>). This package is registered in the
General registry; to install, type `]` in the Julia command prompt, then enter

```julia
pkg> add SDPAFamily
```

Call `SDPAFamily.Optimizer()` to use this wrapper via `MathOptInterface`, which
is an intermediate layer between low-level solvers (such as SDPA-GMP, SDPA-QD,
and SDPA-DD) and high level modelling languages, such as
[JuMP.jl](https://github.com/JuliaOpt/JuMP.jl) and
[Convex.jl](https://github.com/JuliaOpt/Convex.jl/).

Convex.jl 0.13+ supports MathOptInterface and can be used to solve problems with
the solvers from this package.

JuMP currently only supports `Float64` numeric types, which means that problems
can only be specified to 64-bits of precision, and results can only be recovered
at that level of precision, when using JuMP. This is tracked in the issue
[JuMP#2025](https://github.com/JuliaOpt/JuMP.jl/issues/2025).


## Quick Example

Here is a simple optimization problem formulated with Convex.jl:

```julia
using SDPAFamily, LinearAlgebra
using Convex
y = Semidefinite(3)
p = maximize(lambdamin(y), tr(y) <= 5; numeric_type = BigFloat)
solve!(p, () -> SDPAFamily.Optimizer(presolve=true))
@show p.optval
```

See the documentation linked above for troubleshooting help and usage
information.
