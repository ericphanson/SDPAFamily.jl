using SDPAFamily
using Test

# https://github.com/ericphanson/SDPAFamily.jl/pull/36#issuecomment-682117337
using Libdl
libs = filter!(lib -> occursin("fortran", lib), dllist())
@info libs

const variants = (:sdpa, :sdpa_dd, :sdpa_qd, :sdpa_gmp)

@testset "SDPAFamily" begin

    @info "Starting testset `General utilities`"
    @testset "General utilities" begin
        include("presolve.jl")
        #include("status_test.jl") # FIXME
        include("variant_test.jl")
        include("attributes.jl")
    end

    include("Convex.jl")

    include("MOI_wrapper.jl")

    @info "Starting testset `High-precision example`"
    @testset "High-precision example" begin
        include("high_precision_test.jl")
    end

    # No solve in the test suite should need the fallback of the presolve.
    @info "Presolve fallbacks in the test suite: $(SDPAFamily.FALLBACK_COUNT[])"
    @test SDPAFamily.FALLBACK_COUNT[] == 0
end
