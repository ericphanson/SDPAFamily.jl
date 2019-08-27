"""
    sdpa_gmp_binary_solve!(m::Optimizer, full_input_path, full_output_path, extra_args::Cmd, redundant_entries)

Calls the binary `sdpa_gmp` to solve SDPA-formatted problem specified in a .dat-s file at `full_input_path`. Results are written into the file at `full_output_path`.
`extra_args` is passed on to the binary as additional options, allowing for e.g. custom parameter files. `redundant_entries` is a sorted vector listing indices of linearly dependent constraint which are already removed by `presolve.jl`. The corresponding decision variables are populated as zeros.

This function returns `m` with solutions already populated from results in the output file.
"""
function sdpa_gmp_binary_solve!(m::Optimizer, full_input_path::String, full_output_path::String; extra_args::Cmd = ``, redundant_entries::Vector = [])
    # extra_args = `-p /home/jiazheng/Downloads/sdpa-gmp-7.1.3/param.sdpa`
    sdpa_gmp_path = m.binary_path
    arg = `-ds $full_input_path -o $full_output_path $extra_args`
    if !m.silent
        run(`$sdpa_gmp_path $arg`)
    else
        out = devnull
        run(pipeline(`$sdpa_gmp_path $arg`, stdout = out))
    end
    read_results!(m, full_output_path, redundant_entries);
    return m
end
