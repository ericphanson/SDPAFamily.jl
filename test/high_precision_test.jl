using Convex
using Test
time = @elapsed using SDPAFamily
@info "`using SDPAFamily`" time

@testset "High precision test" begin

    E12, E21 = ComplexVariable(2, 2), ComplexVariable(2, 2)
    s1, s2 = [big"0.25" big"-0.25"*im; big"0.25"*im big"0.25"], [big"0.5" big"0.0"; big"0.0" big"0.0"]
    p = Problem{BigFloat}(:minimize, real(tr(E12 * (s1 + 2 * s2) + E21 * (s2 + 2 * s1))), [E12 ⪰ 0, E21 ⪰ 0, E12 + E21 == Diagonal(ones(2)) ])

    time = @elapsed solve!(p, SDPAFamily.Optimizer(silent = true, presolve = true, variant = :sdpa))
    @test p.optval ≈ (6 - sqrt(big"2.0"))/4 atol=1e-6
    error = abs(p.optval - (6 - sqrt(big"2.0"))/4)
    @info "SDPA solved the problem with an absolute error of " error time

    time = @elapsed solve!(p, SDPAFamily.Optimizer(silent = false, presolve = true, variant = :sdpa_dd))
    @test p.optval ≈ (6 - sqrt(big"2.0"))/4 atol=1e-9
    error = abs(p.optval - (6 - sqrt(big"2.0"))/4)
    @info "SDPA-dd solved the problem with an absolute error of " error time

    time = @elapsed solve!(p, SDPAFamily.Optimizer(silent = true, presolve = true, variant = :sdpa_qd))
    error = abs(p.optval - (6 - sqrt(big"2.0"))/4)
    @test p.optval ≈ (6 - sqrt(big"2.0"))/4 atol=1e-15
    @info "SDPA-qd solved the problem with an absolute error of " error time

    time = @elapsed solve!(p, SDPAFamily.Optimizer(silent = true, presolve = true, variant = :sdpa_gmp))
    @test p.optval ≈ (6 - sqrt(big"2.0"))/4 atol=1e-30
    error = abs(p.optval - (6 - sqrt(big"2.0"))/4)
    @info "SDPA-gmp solved the problem with an absolute error of " error time

end
