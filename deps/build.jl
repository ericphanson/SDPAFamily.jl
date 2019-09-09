include("install_bb_qd.jl")
include("install_bb_sdpa.jl")
include("install_bb_sdpa_high_precision.jl")
include("install_custom_WSL_binaries.jl")

using BinaryProvider # requires BinaryProvider 0.3.0 or later

# Parse some basic command-line arguments
const verbose = "--verbose" in ARGS
const prefix = Prefix(get([a for a in ARGS if a != "--verbose"], 1, joinpath(@__DIR__, "usr")))

## Check for WSL or custom library

custom_library = false
HAS_WSL = false

# if the user sets the environmental variable, we take that as the path
# and don't do Product downloads or OS checks.
if haskey(ENV,"JULIA_SDPA_GMP_PATH")
    custom_library = true
    if haskey(ENV,"JULIA_SDPA_GMP_WSL")
        if uppercase(ENV["JULIA_SDPA_GMP_WSL"]) == "TRUE"
            HAS_WSL = true
        end
    end
end

if !custom_library
    # try to detect if we are *within* WSL
    # https://github.com/microsoft/WSL/issues/423#issuecomment-221627364
    if Sys.islinux() && occursin("WSL", read(`cat /proc/version`, String))
        install_wsl_binary = true
    else
        install_wsl_binary = false
    end

    if Sys.iswindows()
        # try to see is `sdpa_gmp` is installed on WSL
        if !isempty(Sys.which("wsl"))
            @info "Windows subsystem for Linux detected. Using WSL-compiled binary."
            install_wsl_binary = true
            HAS_WSL = true

        else
            @info "SDPA-GMP does not directly support Windows, and requires Windows Subsystem for Linux (WSL). "
            error("WSL was not detected.")
        end
    end
end

if custom_library
    products = Product[ExecutableProduct(ENV["JULIA_SDPA_GMP_PATH"], "sdpa_gmp", :sdpa_gmp) ]
elseif install_wsl_binary
    products = install_custom_WSL_binaries(prefix, verbose)
else
    products = install_bb_qd(prefix, verbose)
    append!(products, install_bb_sdpa_high_precision(prefix, verbose))
    append!(products, install_bb_sdpa(prefix, verbose)) # install plain SDPA binary too, for comparisons
end

# Write out a deps.jl file that will contain mappings for our products
deps_file_path = joinpath(@__DIR__, "deps.jl")

write_deps_file(deps_file_path, products, verbose = verbose)

# add WSL info
open(deps_file_path, "a") do io
    write(io, """

    const HAS_WSL = $HAS_WSL
    if $install_wsl_binary || $custom_library
        const libsdpa = ""
    end
    """)
end
