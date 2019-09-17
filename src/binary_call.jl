"""
    sdpa_gmp_binary_solve!(m::Optimizer, full_input_path, full_output_path; redundant_entries)

Calls the binary `sdpa_gmp` to solve SDPA-formatted problem specified in a .dat-s file at `full_input_path`. Results are written into the file at `full_output_path`.
`redundant_entries` is a sorted vector listing indices of linearly dependent constraint which are already removed by `presolve.jl`. The corresponding decision variables are populated as zeros.

This function returns `m` with solutions already populated from results in the output file.
"""
function sdpa_gmp_binary_solve!(m::Optimizer, full_input_path::String, full_output_path::String; redundant_entries::Vector = [])
    read_path = full_output_path
    if m.use_WSL
        full_input_path = WSLize_path(full_input_path)
        full_output_path = WSLize_path(full_output_path)
        if !m.verbosity != SILENT
            @info "Redirecting to WSL environment."
        end
    end
    if startswith(m.params_path, "-pt ") && length(m.params_path) == 5
        arg = `-ds $full_input_path -o $full_output_path $(split(m.params_path))`
    else
        arg = `-ds $full_input_path -o $full_output_path -p $(m.params_path)`
    end
    if m.use_WSL
        wsl_binary_path = dirname(normpath(m.binary_path))
        error_messages, miss = cd(wsl_binary_path) do
            var = string(m.variant)
            run_binary(`wsl ./$var $arg`, m.verbosity)
        end
    else
        error_messages, miss = withenv([prefix]) do
            run_binary(`$(m.binary_path) $arg`, m.verbosity)
        end
    end

    error_log_path = joinpath(m.tempdir, "errors.log")
    open(error_log_path, "w") do io
        print(io, error_messages)
    end

    if m.verbosity != SILENT
        if m.verbosity == VERBOSE && error_messages != ""
            println("error log: $error_log_path") 
        end

        if miss
            @warn("'cholesky miss condition' warning detected; results may be unreliable. Try `presolve=true`, or see troubleshooting guide.")
        end
    end

    read_results!(m, read_path, redundant_entries);
end


function run_binary(cmd::Cmd, verbosity)
    if verbosity == SILENT
        error_messages, miss = run_and_parse_output(devnull, devnull, cmd)
    elseif verbosity == WARN
        error_messages, miss = run_and_parse_output(stdout, stderr, cmd)
    else
        error_messages, miss = run_and_parse_output(stdout, stderr, cmd; echo = true)
    end

    return error_messages, miss
end


function run_and_parse_output(out_io, err_io, cmd; echo = false)
    buffer = IOBuffer()
    miss = false
    function print_stream(out)
        for line in eachline(out)
            if occursin(" :: ", line)
                if occursin("cholesky miss condition", line)
                    miss=true
                end
                println(out_io, "Warning: $line")
                println(buffer, "Warning: $line")
            elseif echo
                println(out_io, line)
            end
        end
    end

    function error_stream(out)
        for line in eachline(out)
            println(err_io, "error: $line")
            write(buffer, "error: $line")
        end
    end

    out = Pipe()
    err = Pipe()
    process = run(pipeline(cmd, stdout=out, stderr=err), wait=false)
    close(out.in)
    close(err.in)


    s1 = @async print_stream(out)
    s2 = @async error_stream(err)
    wait.((process, s1, s2))
    return String(take!(buffer)), miss
end
