export sdpa_gmp_binary

# requires directory with `sdpa_gmp` binary in your PATH variable
const sdpa_gmp_path = "sdpa_gmp"

"""
    sdpa_gmp_binary(arg)
Execute the given command literal as an argument to sdpa_gmp.

"""

# function  sdpa_gmp_binary(arg::Cmd)
#
#     run(pipeline(`$sdpa_gmp_path $arg`, stdout = devnull))
# end

sdpa_gmp_binary(arg::Cmd) = run(`$sdpa_gmp_path $arg`)
