using Convex
using Test
using SDPA_GMP

@testset "High precision test" begin
    E12, E21 = ComplexVariable(2, 2), ComplexVariable(2, 2)
    s1, s2 = [big"0.25" big"-0.25"*im; big"0.25"*im big"0.25"], [big"0.5" big"0.0"; big"0.0" big"0.0"]
    p = Problem{BigFloat}(:minimize, real(tr(E12 * (s1 + 2 * s2) + E21 * (s2 + 2 * s1))), [E12 ⪰ 0, E21 ⪰ 0, E12 + E21 == Diagonal(ones(2)) ])
    p2 = Problem{Float64}(:minimize, real(tr(E12 * (s1 + 2 * s2) + E21 * (s2 + 2 * s1))), [E12 ⪰ 0, E21 ⪰ 0, E12 + E21 == Diagonal(ones(2)) ])

    solve!(p, SDPA_GMP.Optimizer(silent = false, presolve = true, variant = :dd))
    @test p.optval ≈ (6 - sqrt(big"2.0"))/4 atol=1e-7

    solve!(p, SDPA_GMP.Optimizer(silent = false, presolve = true, variant = :qd))
    @test p.optval ≈ (6 - sqrt(big"2.0"))/4 atol=1e-10

    solve!(p, SDPA_GMP.Optimizer(silent = false, presolve = true, variant = :gmp))
    @test p.optval ≈ (6 - sqrt(big"2.0"))/4 atol=1e-30


end