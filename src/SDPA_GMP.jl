module SDPA_GMP
using LinearAlgebra # for diag()
using DelimitedFiles # for writing SDPA input file
using MathOptInterface
const MOI = MathOptInterface
const MOIB = MOI.Bridges
using BinaryProvider


# The `deps.jl` file defines `HAS_WSL::Bool` and `sdpa_gmp::String`.
#
# `HAS_WSL == true` means we default to `use_WSL = true`
# in SDPA_GMP.Optimizer; otherwise we default to `use_WSL = false`.
# if `use_WSL` is set to true for a given SDPA_GMP.Optimizer, that means
# we default to using `default_params_path_wsl` as the params, and
# in the binary call, we turn the paths into WSL paths.
#
# `sdpa_gmp` is the default path to the binary.
include(joinpath(@__DIR__, "..", "deps", "deps.jl"))

function __init__()
    check_deps()
end

"""
    BB_PATHS::Dict{Symbol,String}

Holds the binary-builder-built paths to the executables for `sdpa_gmp`, `sdpa_dd`, and `sdpa_qd`.
"""
const BB_PATHS = Dict(:gmp => sdpa_gmp, :dd => sdpa_dd, :qd => sdpa_qd, :plain => sdpa)

"""
    WSLize_path(path::String) -> String

This function converts Windows paths for use via WSL.
"""
WSLize_path(path) = replace(path, ":" => "") |> x -> replace(x, "\\" => "/") |> x -> "/mnt/"*x |> lowercase

const default_gmp_params_path = normpath(joinpath(@__DIR__, "..", "deps", "param_gmp.sdpa"))
const default_gmp_params_path_wsl = WSLize_path(default_gmp_params_path)

const default_dd_params_path = normpath(joinpath(@__DIR__, "..", "deps", "param_dd.sdpa"))
const default_dd_params_path_wsl = WSLize_path(default_dd_params_path)

const default_qd_params_path = normpath(joinpath(@__DIR__, "..", "deps", "param_qd.sdpa"))
const default_qd_params_path_wsl = WSLize_path(default_qd_params_path)

const default_plain_params_path = normpath(joinpath(@__DIR__, "..", "deps", "param_plain.sdpa"))
const default_plain_params_path_wsl = WSLize_path(default_plain_params_path)

include("MOI_wrapper.jl")
include("file_io.jl")
include("presolve.jl")
include("binary_call.jl")

end # module
