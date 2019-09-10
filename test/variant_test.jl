using Convex, LinearAlgebra, SDPAFamily


@testset "Simple test with variant $var" for var in (:plain, :gmp, :qd, :dd)
    solver = SDPAFamily.Optimizer{BigFloat}(variant=var, silent=true)
    y = Variable((3, 3), :Semidefinite)
    p = Problem{BigFloat}(:minimize, tr(y), y[2,1]<=4, y[2,2]>=3)
    solve!(p, solver)
    @test p.optval ≈ big"3.0" atol=1e-5
end
