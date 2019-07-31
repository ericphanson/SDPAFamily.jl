module SDPA_GMP
using Convex

export sdpa_gmp_binary_solve


include("write_problem.jl")
include("binary_call.jl")
include("read_results.jl")

function sdpa_gmp_binary_solve(p::Problem{T}, extra_args::Cmd = `-pt 1`) where {T}
    mock = SDOI.MockSDOptimizer{T}()
    temp = mktempdir()
    inputname = "input.dat-s"
    outputname = "output.dat"
    full_input_path = joinpath(temp, inputname)
    full_output_path = joinpath(temp, outputname)
    write_problem(mock, p, full_input_path)
    #
    sdpa_gmp_binary(`$full_input_path $full_output_path $extra_args`)
    # println(read(full_output_path, String))
    read_results(mock, full_output_path)
    return mock
end

end # module
