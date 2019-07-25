module SDPA_GMP

export sdpa_gmp_binary

using BinaryProvider

const sdpa_gmp_path = joinpath(@__DIR__, "..", "deps", "usr", "bin", "sdpa_gmp")

# Load in `deps.jl`, complaining if it does not exist
const depsjl_path = joinpath(@__DIR__, "..", "deps", "deps.jl")
if !isfile(depsjl_path)
    println("Deps path: $depsjl_path")
    error("SDPA_GMP not installed properly, run `] build SDPA_GMP`, restart Julia and try again")
end

include(depsjl_path)

"""
    sdpa_gmp_binary(arg)
Execute the given command literal as an argument to sdpa_gmp.

"""
sdpa_gmp_binary(arg::Cmd) = run(`$sdpa_gmp_path $arg`)

end # module
