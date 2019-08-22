using Convex
using SDPA_GMP
using MathOptInterface
const MOI = MathOptInterface
const MOIU = MOI.Utilities
const MOIB = MOI.Bridges
using SemidefiniteOptInterface
const SDOI = SemidefiniteOptInterface

using SparseArrays
# using SDPA

using LinearAlgebra
using GenericLinearAlgebra

eye(n) = Matrix(big"1.0"*I, n, n)

using SCS

# # x = [1, 2, 3]
# # P = Variable(3, 3)
# # p1 = Problem{BigFloat}(:minimize, matrixfrac(x, P), [P <= big(2)*eye(3), P >= big(0.5) * eye(3)])

# x = Semidefinite(3)
# p = Problem{BigFloat}(:minimize, sumlargesteigs(x, 2), [x >= big(1)])
# #
# # mock1 = SDPA_GMP.sdpa_gmp_binary_solve(p1)
# #
# Y = Variable(5,5)
# X = -1*rand(BigFloat, 5, 5)
# p2 = Problem{BigFloat}(:minimize, tr(Y), [ diag(Y)[2:5] == diag(X)[2:5], Y[1,1] == big(0.0) ])
# #

# x = Variable(1)
# y = Variable(1)
# p2 = Problem{BigFloat}(:maximize, x + y, [ x <= big"1.0" , y <= big"2.0"])

# mock = Convex.solve!(p2, SDPA_GMP.SDPAGMPoptimizer(BigFloat, verbose = true));

# # mock2 = SDPA_GMP.sdpa_gmp_binary_solve(p2)

E12, E21 = ComplexVariable(2, 2), ComplexVariable(2, 2);
s1, s2 = [big"0.25" big"-0.25"*im; big"0.25"*im big"0.25"], [big"0.5" big"0.0"; big"0.0" big"0.0"]

p3 = Problem{BigFloat}(:minimize, tr(real(E12 * (s1 + big"2.0" * s2) + E21 * (s2 + 2 * s1))), [E12 ⪰ 0, E12 + E21 == Diagonal([big"1.0", big"1.0"])])
# #
# # mock3 = SDPA_GMP.sdpa_gmp_binary_solve(p3)
# x = Variable(Positive())
# y = Semidefinite(3)
# p = Problem{BigFloat}(:minimize, nuclearnorm(y), y[2,1]<=4, y[2,2]>=3, y[3,3]<=2)
# solve!(p, SDPA_GMP.SDPAGMPoptimizer(BigFloat, verbose = true));

# x = Variable(4)
# p = Problem{BigFloat}(:minimize, sum(Diagonal(x)), x == [1, 2, 3, 4])
# @test vexity(p) == AffineVexity()
# solve!(p, solver)
# @test p.optval ≈ 10 atol=TOL
# @test all(abs.(evaluate(Diagonal(x)) - Diagonal([1, 2, 3, 4])) .<= TOL)

# x = Variable(4,1)
# p = Problem{BigFloat}(:minimize, dotsort(x, [1,2,3,4]), sum(x) >= 7, x >=0, x<=2, x[4]<=1)

# x = Variable(3)
# p = Problem{BigFloat}(:minimize, norm_1(x), [-2 <= x, x <= 1])

# x = Variable(4, 4)
# p = Problem{BigFloat}(:minimize, sumlargest(x, 2), sumsmallest(x, 4) >= 1)

# A = Semidefinite(2)
# B = [1 0; 0 0]
# ρ = kron(B, A)
# constraints = [partialtrace(ρ, 1, [2; 2]) == BigFloat[0.09942819 0.29923607; 0.29923607 0.90057181], ρ in :SDP]
# p = Problem{BigFloat}(:minimize, Constant(0), constraints)

a = big"2.214123515"+big"0.8143999"*im
x = ComplexVariable()
objective = norm2(a-x)
c1 = real(x)>=0
p = Problem{BigFloat}(:minimize, objective,c1)
mock = solve!(p3, SDPA_GMP.Optimizer{BigFloat}());
solve!(p, ProxSDP.Optimizer(log_verbose = true));
# d = solve!(p, e);
# g = SDOI.mockSDoptimizer(Float64);
# f = SCS.Optimizer();

# g = [1.2; 4.3];
# u = Variable(2);
# objective = norm2(u-g);
# c2 = [];
# r = minimize(objective);



n = 70
M = rand(n,n)+ im*rand(n,n)
M = M + M'
 # now M is hermitian
x = ComplexVariable(n,n)
objective = norm2(M - x)
c1 = x in :SDP
p = Problem{Float64}(:minimize, objective, c1)
mmo = solve!(p3, SDPA_GMP.SDPAGMPoptimizer(BigFloat, verbose = true, normal_sdpa=false));
c, A, b, cones, var_to_ranges, vartypes, conic_constraints = Convex.conic_problem(p3);
