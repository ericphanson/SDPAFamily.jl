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

    time = @elapsed solve!(p, SDPAFamily.Optimizer(silent = true, presolve = true, variant = :sdpa_dd))
    @test p.optval ≈ (6 - sqrt(big"2.0"))/4 atol=1e-15
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


# This example is in the documentation, so let's make sure it continues to work.
@testset "Very high precision test" begin
    opt = SDPAFamily.Optimizer(
        presolve = true,
        params = (  epsilonStar = 1e-200, # constraint tolerance
                    epsilonDash = 1e-200, # normalized duality gap tolerance
                    precision = 2000 # arithmetric precision used in sdpa_gmp
                ))

    setprecision(2000) do
        ρ₁ = Complex{BigFloat}[1 0; 0 0]
        ρ₂ = (1//2)*Complex{BigFloat}[1 -im; im 1]
        E₁ = ComplexVariable(2, 2);
        E₂ = ComplexVariable(2, 2);
        problem = maximize( real((1//2)*tr(ρ₁*E₁) + (1//2)*tr(ρ₂*E₂)),
                            [E₁ ⪰ 0, E₂ ⪰ 0, E₁ + E₂ == Diagonal(ones(2))];
                            numeric_type = BigFloat );
        solve!(problem, opt)
        p_guess = 1//2 + 1/(2*sqrt(big(2)))
        @test problem.optval ≈ p_guess atol=1e-200
    end
end
