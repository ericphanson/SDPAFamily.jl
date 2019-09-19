module SDPAFamily
using LinearAlgebra # for diag()
using DelimitedFiles # for writing SDPA input file
using MathOptInterface
const MOI = MathOptInterface
const MOIB = MOI.Bridges
using BinaryProvider

"""
Possible verbosity levels of an `SDPAFamily.Optimizer`.

Options are `SILENT`, `WARN`, or `VERBOSE`.
"""
@enum Verbosity SILENT WARN VERBOSE

"""
Possible the parameter settings of an `SDPAFamily.Optimizer`.
One can also pass a path to the `params` keyword argument to use
a custom parameter file.

Options are `DEFAULT`, `UNSTABLE_BUT_FAST`, or `STABLE_BUT_SLOW`.
"""
@enum ParamsSetting DEFAULT UNSTABLE_BUT_FAST STABLE_BUT_SLOW


# The `deps.jl` file defines `HAS_WSL::Bool` and `sdpa_gmp::String`.
#
# `HAS_WSL == true` means we default to `use_WSL = true`
# in SDPAFamily.Optimizer; otherwise we default to `use_WSL = false`.
# if `use_WSL` is set to true for a given SDPAFamily.Optimizer, that means
# we default to using `default_params_path_wsl` as the params, and
# in the binary call, we turn the paths into WSL paths.
#
# `sdpa_gmp` is the default path to the binary.
include(joinpath(@__DIR__, "..", "deps", "deps.jl"))

const prefix = Prefix(joinpath(@__DIR__, "..", "deps", "usr"))

function __init__()
    check_deps()
end

"""
    BB_PATHS::Dict{Symbol,String}

Holds the binary-builder-built paths to the executables for `sdpa_gmp`, `sdpa_dd`, and `sdpa_qd`.
"""
const BB_PATHS = Dict(:sdpa_gmp => sdpa_gmp, :sdpa_dd => sdpa_dd, :sdpa_qd => sdpa_qd, :sdpa => sdpa)

const default_params_path = Dict(
    :sdpa_gmp =>  normpath(joinpath(@__DIR__, "..", "deps", "param_gmp.sdpa")),
    :sdpa_dd => normpath(joinpath(@__DIR__, "..", "deps", "param_dd.sdpa")),
    :sdpa_qd => normpath(joinpath(@__DIR__, "..", "deps", "param_qd.sdpa")),
    :sdpa => normpath(joinpath(@__DIR__, "..", "deps", "param_plain.sdpa")),
    :sdpa_gmp_float64 => normpath(joinpath(@__DIR__, "..", "deps", "param_gmp_float64.sdpa"))
)

"""
    WSLize_path(path::String) -> String

This function converts Windows paths for use via WSL.
"""
WSLize_path(path) = replace(path, ":" => "") |> x -> replace(x, "\\" => "/") |> x -> "/mnt/"*x |> lowercase

include("params.jl")
include("MOI_wrapper.jl")
include("file_io.jl")
include("presolve.jl")
include("binary_call.jl")

end # module
