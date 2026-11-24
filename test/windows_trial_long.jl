# The package's MOI and Convex tests for `:sdpa_gmp`, run by `.github/workflows/windows_trial.yml`
# after `windows_trial.jl` passes.
using Test

const variants = (:sdpa_gmp,)

@testset "Windows trial, long tests" begin
    include("MOI_wrapper.jl")
    include("Convex.jl")
end
