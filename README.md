# SDPAFamily

[![CI](https://github.com/ericphanson/SDPAFamily.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/ericphanson/SDPAFamily.jl/actions/workflows/ci.yml)
[![Codecov](https://codecov.io/gh/ericphanson/SDPAFamily.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/ericphanson/SDPAFamily.jl)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://ericphanson.github.io/SDPAFamily.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://ericphanson.github.io/SDPAFamily.jl/dev)

Note: in v0.3.0, which updated the package for MOI v1, the `presolve` option
had no effect, even when set to `true`. Without a presolve, SDPA can stop at the
iteration limit or return inaccurate results for problems with equality
constraints or linearly dependent constraints. The presolve is restored, and
enabled by default, in the development version; see issue
[#67](https://github.com/ericphanson/SDPAFamily.jl/issues/67) and the
[documentation](https://ericphanson.github.io/SDPAFamily.jl/dev/usage/#Presolve).

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
p = maximize(eigmin(y), tr(y) <= 5; numeric_type = BigFloat)
solve!(p, () -> SDPAFamily.Optimizer())
@show p.optval
```

See the documentation linked above for troubleshooting help and usage
information.
