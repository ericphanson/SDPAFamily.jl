# SDPAFamily

This package provides an interface to using SDPA-GMP, SDPA-DD, and SDPA-QD in
Julia (<http://sdpa.sourceforge.net>). Call `SDPAFamily.Optimizer()` to use this
wrapper via `MathOptInterface`, which is an intermediate layer between low-level
solvers (such as SDPA-GMP, SDPA-QD, and SDPA-DD) and high level modelling
languages, such as [JuMP.jl](https://github.com/JuliaOpt/JuMP.jl) and
[Convex.jl](https://github.com/JuliaOpt/Convex.jl/).

To use the usual (non high-precision) SDPA solver, see
[SDPA.jl](https://github.com/JuliaOpt/SDPA.jl).


Convex.jl 0.13+ supports MathOptInterface and can be used to solve problems with
the solvers from this package.

JuMP currently only supports `Float64` numeric types, which means that problems
can only be specified to 64-bits of precision, and results can only be recovered
at that level of precision, when using JuMP. This is tracked in the issue
[JuMP#2025](https://github.com/JuliaOpt/JuMP.jl/issues/2025).


## Credits

This package was developed by Jiazheng Zhu in the summer of 2019, under the
supervision of Eric Hanson.
