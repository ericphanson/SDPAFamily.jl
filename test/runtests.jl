using SDPA_GMP
using Test

@testset "SDPA_GMP.jl" begin
    sdpa_gmp_binary(`-o example1.result -dd example1.dat-s`)
end
