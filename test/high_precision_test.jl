using Convex
using Test
time = @elapsed using SDPA_GMP
@info "`using SDPA_GMP`" time

@testset "High precision test" begin
    E12, E21 = ComplexVariable(2, 2), ComplexVariable(2, 2)
    s1, s2 = [big"0.25" big"-0.25"*im; big"0.25"*im big"0.25"], [big"0.5" big"0.0"; big"0.0" big"0.0"]
    p = Problem{BigFloat}(:minimize, real(tr(E12 * (s1 + 2 * s2) + E21 * (s2 + 2 * s1))), [E12 ⪰ 0, E21 ⪰ 0, E12 + E21 == Diagonal(ones(2)) ])

    time = @elapsed solve!(p, SDPA_GMP.Optimizer(silent = true, presolve = true, variant = :plain))
    @test p.optval ≈ (6 - sqrt(big"2.0"))/4 atol=1e-10
    error = abs(p.optval - (6 - sqrt(big"2.0"))/4)
    @info "SDPA solved the problem with an absolute error of " error time

    time = @elapsed solve!(p, SDPA_GMP.Optimizer(silent = true, presolve = true, variant = :dd))
    @test p.optval ≈ (6 - sqrt(big"2.0"))/4 atol=1e-10
    error = abs(p.optval - (6 - sqrt(big"2.0"))/4)
    @info "SDPA-dd solved the problem with an absolute error of " error time

    time = @elapsed solve!(p, SDPA_GMP.Optimizer(silent = true, presolve = true, variant = :qd))
    error = abs(p.optval - (6 - sqrt(big"2.0"))/4)
    @test p.optval ≈ (6 - sqrt(big"2.0"))/4 atol=1e-10
    @info "SDPA-qd solved the problem with an absolute error of " error time

    time = @elapsed solve!(p, SDPA_GMP.Optimizer(silent = true, presolve = true, variant = :gmp))
    @test p.optval ≈ (6 - sqrt(big"2.0"))/4 atol=1e-30
    error = abs(p.optval - (6 - sqrt(big"2.0"))/4)
    @info "SDPA-gmp solved the problem with an absolute error of " error time


end