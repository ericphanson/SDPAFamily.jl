module SDPA_GMP
using LinearAlgebra # for diag()
using DelimitedFiles # for writing SDPA input file
using MathOptInterface
const MOI = MathOptInterface
const MOIB = MOI.Bridges
using BinaryProvider


# The `deps.jl` file defines `has_WSL::Bool` and `sdpa_gmp::String`.
#
# `has_WSL == true` means we default to `use_WSL = true`
# in SDPA_GMP.Optimizer; otherwise we default to `use_WSL = false`.
# if `use_WSL` is set to true for a given SDPA_GMP.Optimizer, that means
# we default to using `default_params_path_wsl` as the params, and
# in the binary call, we turn the paths into WSL paths.
#
# `sdpa_gmp` is the default path to the binary.
include(joinpath(@__DIR__, "..", "deps", "deps.jl"))

"""
    WSLize_path(path::String) -> String

This function converts Windows paths for use via WSL.
"""
WSLize_path(path) = replace(path, ":" => "") |> x -> replace(x, "\\" => "/") |> x -> "/mnt/"*x |> lowercase

const default_params_path = normpath(joinpath(@__DIR__, "..", "deps", "param.sdpa"))
const default_params_path_wsl = WSLize_path(default_params_path)


include("MOI_wrapper.jl")
include("file_io.jl")
include("presolve.jl")
include("binary_call.jl")

end # module
