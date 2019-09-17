using Test, SDPAFamily
using Convex, LinearAlgebra

params_file = joinpath(@__DIR__, "custom_params.sdpa")
setprecision(800) do
    E12, E21 = ComplexVariable(2, 2), ComplexVariable(2, 2)
    s1, s2 = [big"0.25" big"-0.25"*im; big"0.25"*im big"0.25"], [big"0.5" big"0.0"; big"0.0" big"0.0"]
    p = Problem{BigFloat}(:minimize, real(tr(E12 * (s1 + 2 * s2) + E21 * (s2 + 2 * s1))), [E12 ⪰ 0, E21 ⪰ 0, E12 + E21 == Diagonal(ones(2)) ])
    opt = SDPAFamily.Optimizer(params = params_file, presolve = true, variant = :sdpa_gmp)
    time = @elapsed solve!(p, opt)

    correct_val =  (6 - sqrt(big"2.0"))/4
    @test p.optval ≈ correct_val atol=1e-80
end
