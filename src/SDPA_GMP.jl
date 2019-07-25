module SDPA_GMP

using BinaryProvider

const libpath = joinpath(@__DIR__, "..", "deps", "usr", "lib")

if Sys.iswindows()
    const execenv = ("PATH" => string(libpath, ";", Sys.BINDIR))
elseif Sys.isapple()
    const execenv = ("DYLD_LIBRARY_PATH" => libpath)
else
    const execenv = ("LD_LIBRARY_PATH" => libpath)
end


# Load in `deps.jl`, complaining if it does not exist
const depsjl_path = joinpath(@__DIR__, "..", "deps", "deps.jl")
if !isfile(depsjl_path)
    println("Deps path: $depsjl_path")
    error("FFMPEG not installed properly, run `] build FFMPEG`, restart Julia and try again")
end

include(depsjl_path)

"""
    sdpa_gmp(arg)
Execute the given command literal as an argument to sdpa_gmp.

"""
function sdpa_gmp(arg::Cmd)
    withenv(execenv) do
            Base.run(`$command $arg`)
        end
end
end # module
