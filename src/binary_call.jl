"""
    sdpa_gmp_binary_solve!(m::Optimizer, full_input_path, full_output_path; extra_args::Cmd, redundant_entries)

Calls the binary `sdpa_gmp` to solve SDPA-formatted problem specified in a .dat-s file at `full_input_path`. Results are written into the file at `full_output_path`.
`extra_args` is passed on to the binary as additional options, allowing for e.g. custom parameter files. `redundant_entries` is a sorted vector listing indices of linearly dependent constraint which are already removed by `presolve.jl`. The corresponding decision variables are populated as zeros.

This function returns `m` with solutions already populated from results in the output file.
"""
function sdpa_gmp_binary_solve!(m::Optimizer, full_input_path::String, full_output_path::String; extra_args::Cmd = ``, redundant_entries::Vector = [])
    read_path = full_output_path
    if m.use_WSL
        full_input_path = WSLize_path(full_input_path)
        full_output_path = WSLize_path(full_output_path)
        if !m.silent
            @info "Redirecting to sdpa_gmp in WSL environment."
        end
    end
    arg = `-ds $full_input_path -o $full_output_path -p $(m.params_path) $extra_args`
    if m.use_WSL
        wsl_binary_path = dirname(m.binary_path)
        cd(wsl_binary_path) do
            withenv([prefix]) do
                run(pipeline(`wsl sdpa_gmp $arg`, stdout = m.silent ? devnull : stdout))
            end
        end
    else
        withenv([prefix]) do
            run(pipeline(`$(m.binary_path) $arg`, stdout = m.silent ? devnull : stdout))
        end
    end
    read_results!(m, read_path, redundant_entries);
end
