```@example
using SDPAFamily
using Pkg
Pkg.add(PackageSpec(name="Convex", url="https://github.com/ericphanson/Convex.jl", rev="MathOptInterface"));
using Convex

E12, E21 = ComplexVariable(2, 2), ComplexVariable(2, 2)
s1, s2 = [big"0.25" big"-0.25"*im; big"0.25"*im big"0.25"], [big"0.5" big"0.0"; big"0.0" big"0.0"]
p = Problem{BigFloat}(:minimize, real(tr(E12 * (s1 + 2 * s2) + E21 * (s2 + 2 * s1))), [E12 ⪰ 0, E21 ⪰ 0, E12 + E21 == Diagonal(ones(2)) ])

time = @elapsed solve!(p, SDPAFamily.Optimizer(silent = true, presolve = true, variant = :sdpa))
error = abs(p.optval - (6 - sqrt(big"2.0"))/4)
@info "SDPA solved the problem with an absolute error of " error time

time = @elapsed solve!(p, SDPAFamily.Optimizer(silent = true, presolve = true, variant = :sdpa_dd))
error = abs(p.optval - (6 - sqrt(big"2.0"))/4)
@info "SDPA-dd solved the problem with an absolute error of " error time

time = @elapsed solve!(p, SDPAFamily.Optimizer(silent = true, presolve = true, variant = :sdpa_qd))
error = abs(p.optval - (6 - sqrt(big"2.0"))/4)
@info "SDPA-qd solved the problem with an absolute error of " error time

time = @elapsed solve!(p, SDPAFamily.Optimizer(silent = true, presolve = true, variant = :sdpa_gmp))
error = abs(p.optval - (6 - sqrt(big"2.0"))/4)
@info "SDPA-gmp solved the problem with an absolute error of " error time
```
