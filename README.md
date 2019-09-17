# SDPAFamily

[![Build Status](https://travis-ci.com/ericphanson/SDPAFamily.jl.svg?branch=master)](https://travis-ci.com/ericphanson/SDPAFamily.jl)
[![Codecov](https://codecov.io/gh/ericphanson/SDPAFamily.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/ericphanson/SDPAFamily.jl)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://ericphanson.github.io/SDPAFamily.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://ericphanson.github.io/SDPAFamily.jl/dev)

An interface to using SDPA-GMP, SDPA-DD, and SDPA-QD in Julia (<http://sdpa.sourceforge.net>). Call `SDPAFamily.Optimizer()` to use this wrapper via `MathOptInterface`, which is an intermediate layer between low-level solvers (such as SDPA-GMP, SDPA-QD, and SDPA-DD) and high level modelling languages, such as [JuMP.jl](https://github.com/JuliaOpt/JuMP.jl) and [Convex.jl](https://github.com/JuliaOpt/Convex.jl/).

JuMP currently only supports `Float64` numeric types, which means that problems can only be specified to 64-bits of precision, and results can only be recovered at that level of precision, when using JuMP. This is tracked in the issue [JuMP#2025](https://github.com/JuliaOpt/JuMP.jl/issues/2025).

Convex.jl does not yet officially support MathOptInterface; this issue is tracked at [Convex.jl#262](https://github.com/JuliaOpt/Convex.jl/issues/262). However, there is a work-in-progress branch which can be added to your Julia environment via

```julia
] add https://github.com/ericphanson/Convex.jl#MathOptInterface
```

which can be used to solve problems with the solvers from this package.

## Quick Example

Here is a simple optimization problem formulated with Convex.jl:

```julia
using SDPAFamily, LinearAlgebra
using Convex # ] add https://github.com/ericphanson/Convex.jl#MathOptInterface
y = Semidefinite(3)
p = maximize(lambdamin(y), tr(y) <= 5; numeric_type = BigFloat)
solve!(p, SDPAFamily.Optimizer(presolve=true))
@show p.optval
```

See the documentation linked above for troubleshooting help and usage information.
