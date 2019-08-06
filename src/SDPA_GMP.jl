module SDPA_GMP
using MathOptInterface
MOI = MathOptInterface
using LinearAlgebra # for diag()
MOIU = MOI.Utilities
MOIB = MOI.Bridges
const sdpa_gmp_path = "sdpa_gmp" # make sure sdpa_gmp binary is in the system PATH

export sdpa_gmp_binary_solve!, read_results!

include("moi_wrap.jl")
include("read_results.jl")


"""
    sdpa_gmp_binary(arg; verbose::Bool)
Execute the given command literal as the only argument to sdpa_gmp.

"""
function  sdpa_gmp_binary(arg::Cmd; verbose = false)
    if verbose
        sdpa_gmp_binary(arg::Cmd) = run(`$sdpa_gmp_path $arg`)
    else
        run(pipeline(`$sdpa_gmp_path $arg`, stdout = devnull))
    end
end

"""
    sdpa_gmp_binary_solve!(m::SDPAGMPOptimizer, full_input_path, full_output_path, extra_args::Cmd)

Calls the binary `sdpa_gmp` to solve SDPA-formatted problem specified in a .dat-s file at `full_input_path`. Results are written into the file at `full_output_path`.
`extra_args` is passed on to the binary as additional options, allowing for e.g. custom parameter files.

This function returns `m` with solutions already populated from results in the output file.
"""
function sdpa_gmp_binary_solve!(m::SDPAGMPOptimizer, full_input_path::String, full_output_path::String, extra_args::Cmd = `-p /home/jiazheng/Downloads/sdpa-gmp-7.1.3/param.sdpa`)
    sdpa_gmp_binary(`-ds $full_input_path -o $full_output_path $extra_args`);
    # println(read(full_output_path, String))
    read_results!(m, full_output_path);
    return m
end

end # module
