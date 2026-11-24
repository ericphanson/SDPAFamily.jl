# Checks for the trial Windows build of SDPA_GMP_jll; `.github/workflows/windows_trial.yml` runs them.
using Test, Downloads, LinearAlgebra
using Convex, SDPAFamily
using SDPAFamily: SDPA_GMP_jll

@testset "Windows trial" begin
    # Upstream's example1.dat-s has LF line endings; also try it with CRLF.
    @testset "sdpa_gmp.exe solves example1.dat-s with $eol line endings" for eol in ("LF", "CRLF")
        dir = mktempdir()
        input = Downloads.download("https://raw.githubusercontent.com/nakatamaho/sdpa-gmp/ca110db5ea1cc46e811b70dfea9cbb25db74448d/example1.dat-s",
                                   joinpath(dir, "example1.dat-s"))
        eol == "CRLF" && write(input, replace(read(input, String), "\n" => "\r\n"))
        output = joinpath(dir, "example1.out")
        base = SDPA_GMP_jll.sdpa_gmp()
        run(Cmd(`$base -ds $input -o $output -pt 0`; env = base.env))
        lines = readlines(output)
        @test occursin("pdOPT", only(filter(startswith("phase.value"), lines)))
        objval = parse(BigFloat, split(only(filter(startswith("objValPrimal"), lines)))[3])
        @test objval ≈ -41.9 atol = 1e-10
    end

    @testset "README example" begin
        y = Semidefinite(3)
        p = maximize(eigmin(y), tr(y) <= 5; numeric_type = BigFloat)
        solve!(p, () -> SDPAFamily.Optimizer(presolve = true, silent = true))
        @test p.optval ≈ big"5" / 3 atol = 1e-10
    end

    @testset "variant_test.jl problem" begin
        y = Variable((3, 3), :Semidefinite)
        p = Problem{BigFloat}(:minimize, tr(y), y[2, 1] <= 4, y[2, 2] >= 3)
        solve!(p, () -> SDPAFamily.Optimizer{BigFloat}(variant = :sdpa_gmp, silent = true))
        @test p.optval ≈ big"3.0" atol = 1e-5
    end
end
