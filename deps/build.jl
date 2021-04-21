include("install_bb_qd.jl")
include("install_bb_sdpa_high_precision.jl")
include("install_custom_WSL_binaries.jl")

using BinaryProvider # requires BinaryProvider 0.3.0 or later
using Scratch, UUIDs

const SDPAFamily_UUID = UUID("bfe18334-aefd-11e9-1109-4bf2b15a5b91")

# Key off of Julia version to avoid issues like <https://github.com/ericphanson/SDPAFamily.jl/issues/29#issue-549039097>
# why key off `Sys.KERNEL`? I once had some strange setup with WSL2 where I had windows
# and linux runtimes looking at the same filesystem. This might help...
const DEPS_DIR = get_scratch!(SDPAFamily_UUID, string("build_julia_", VERSION.major, ".", VERSION.minor, "_", Sys.KERNEL))

# Parse some basic command-line arguments
const verbose = "--verbose" in ARGS
const prefix = Prefix(get([a for a in ARGS if a != "--verbose"], 1, joinpath(DEPS_DIR, "usr")))

## Check for WSL or custom library

custom_library = Dict(:sdpa_gmp => false, :sdpa_qd => false, :sdpa_dd => false, :sdpa => false)
HAS_WSL = Dict(:sdpa_gmp => false, :sdpa_qd => false, :sdpa_dd => false, :sdpa => false)
path_names = Dict(:sdpa_gmp => "JULIA_SDPA_GMP_PATH", :sdpa_dd => "JULIA_SDPA_DD_PATH", :sdpa_qd => "JULIA_SDPA_QD_PATH", :sdpa => "JULIA_SDPA_PLAIN_PATH")

# if the user sets the environmental variable, we take that as the path
# and don't do Product downloads or OS checks.
if haskey(ENV,"JULIA_SDPA_GMP_PATH")
    custom_library[:sdpa_gmp] = true
    if haskey(ENV,"JULIA_SDPA_GMP_WSL")
        if uppercase(ENV["JULIA_SDPA_GMP_WSL"]) == "TRUE"
            HAS_WSL[:sdpa_gmp] = true
        end
    end
end
if haskey(ENV,"JULIA_SDPA_DD_PATH")
    custom_library[:sdpa_dd] = true
    if haskey(ENV,"JULIA_SDPA_DD_WSL")
        if uppercase(ENV["JULIA_SDPA_DD_WSL"]) == "TRUE"
            HAS_WSL[:sdpa_dd] = true
        end
    end
end
if haskey(ENV,"JULIA_SDPA_QD_PATH")
    custom_library[:sdpa_qd] = true
    if haskey(ENV,"JULIA_SDPA_QD_WSL")
        if uppercase(ENV["JULIA_SDPA_QD_WSL"]) == "TRUE"
            HAS_WSL[:sdpa_qd] = true
        end
    end
end
if haskey(ENV,"JULIA_SDPA_PLAIN_PATH")
    custom_library[:sdpa] = true
    if haskey(ENV,"JULIA_SDPA_PLAIN_WSL")
        if uppercase(ENV["JULIA_SDPA_PLAIN_WSL"]) == "TRUE"
            HAS_WSL[:sdpa] = true
        end
    end
end

products = Product[]

for var in [:sdpa_gmp, :sdpa_qd, :sdpa_dd]
    if custom_library[var]
        push!(products, FileProduct(Prefix(ENV[path_names[var]]), string(var), var))
    else
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
                HAS_WSL[var] = true

            else
                @info "SDPA-GMP does not directly support Windows, and requires Windows Subsystem for Linux (WSL). "
                error("WSL was not detected.")
            end
        end
        if install_wsl_binary
            append!(products, install_custom_WSL_binaries(prefix, verbose, [var]))
        else
            append!(products, install_bb_sdpa_high_precision(prefix, verbose, [var]))
            append!(products, install_bb_qd(prefix, verbose))
        end
    end
end

# Write out a deps.jl file that will contain mappings for our products
deps_file_path = joinpath(DEPS_DIR, "deps.jl")

write_deps_file(deps_file_path, products, verbose = verbose)

# add WSL info
open(deps_file_path, "a") do io
    write(io, """

    const HAS_WSL = Dict{Symbol,Bool}($HAS_WSL)

    """)
end
