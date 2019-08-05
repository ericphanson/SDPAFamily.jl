module SDPA_GMP
using Convex
using MathOptInterface
MOI = MathOptInterface
using SemidefiniteOptInterface
SDOI = SemidefiniteOptInterface
export sdpa_gmp_binary_solve

include("moi_wrap.jl")
include("write_problem.jl")
include("binary_call.jl");
include("read_results.jl")

function sdpa_gmp_binary_solve(m::SDPAGMPOptimizer, full_input_path::String, full_output_path::String, extra_args::Cmd = `-p /home/jiazheng/Downloads/sdpa-gmp-7.1.3/param.sdpa`)
    # mock = SDOI.MockSDOptimizer{T}()
    # temp = mktempdir()
    # inputname = "input.dat-s"
    # outputname = "output.dat"
    # full_input_path = joinpath(temp, inputname)
    # full_output_path = joinpath(temp, outputname)
    # write_problem(mock, p, full_input_path)
    #
    sdpa_gmp_binary(`-ds $full_input_path -o $full_output_path $extra_args`);
    # println(read(full_output_path, String))
    read_results(m, full_output_path);
    return m
end

end # module
