module SDPA_GMP

const sdpa_gmp_path = joinpath("sdpa_gmp")


"""
    sdpa_gmp_binary(arg)
Execute the given command literal as an argument to sdpa_gmp.

"""
function sdpa_gmp_binary(arg::Cmd)
    Base.run(`$sdpa_gmp_path $arg`) # needs `sdpa_gmp` in your PATH
end


end # module
