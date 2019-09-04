"""
    sdpa_gmp_binary_solve!(m::Optimizer, full_input_path, full_output_path; extra_args::Cmd, redundant_entries)

Calls the binary `sdpa_gmp` to solve SDPA-formatted problem specified in a .dat-s file at `full_input_path`. Results are written into the file at `full_output_path`.
`extra_args` is passed on to the binary as additional options, allowing for e.g. custom parameter files. `redundant_entries` is a sorted vector listing indices of linearly dependent constraint which are already removed by `presolve.jl`. The corresponding decision variables are populated as zeros.

This function returns `m` with solutions already populated from results in the output file.
"""
function sdpa_gmp_binary_solve!(m::Optimizer, full_input_path::String, full_output_path::String; extra_args::Cmd = ``, redundant_entries::Vector = [])
    if m.use_WSL
        full_input_path = replace(full_input_path, ":" => "") |> x -> replace(x, "\\" => "/") |> x -> "/mnt/"*x |> x -> lowercase(x)
        full_output_path = replace(full_output_path, ":" => "") |> x -> replace(x, "\\" => "/") |> x -> "/mnt/"*x |> x -> lowercase(x)
        if !m.silent
            @info "Redirecting to sdpa_gmp in WSL environment."
        end
    end
    arg = `-ds $full_input_path -o $full_output_path -p $(m.params_path) $extra_args`
    if m.use_WSL
        wsl_binary_path = replace(m.binary_path, ":" => "") |> x -> replace(x, "\\" => "/") |> x -> "/mnt/"*x |> x -> lowercase(x)
        run(pipeline(`wsl $wsl_binary_path $arg`, stdout = m.silent ? devnull : stdout))
    else
        run(pipeline(`$(m.binary_path) $arg`, stdout = m.silent ? devnull : stdout))
    end
    read_results!(m, full_output_path, redundant_entries);
end
