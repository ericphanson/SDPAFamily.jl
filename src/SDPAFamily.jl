module SDPAFamily
using SDPA_jll: SDPA_jll
using SDPA_GMP_jll: SDPA_GMP_jll
using SDPA_DD_jll: SDPA_DD_jll
using SDPA_QD_jll: SDPA_QD_jll
using LinearAlgebra # for diag()
using DelimitedFiles # for writing SDPA input file
import MathOptInterface as MOI
using Scratch

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

"""
    default_binary_path(variant::Symbol) -> String

Path to a user-supplied executable for `variant`, or `""` to use the executable from the
corresponding `_jll` package. A custom executable is selected by setting the environment
variable `JULIA_SDPA_GMP_PATH`, `JULIA_SDPA_DD_PATH`, `JULIA_SDPA_QD_PATH` or
`JULIA_SDPA_PLAIN_PATH` to the directory that contains it.
"""
function default_binary_path(variant::Symbol)
    var = PATH_ENV_VARS[variant]
    haskey(ENV, var) || return ""
    return joinpath(ENV[var], string(variant))
end

const PATH_ENV_VARS = Dict(:sdpa_gmp => "JULIA_SDPA_GMP_PATH", :sdpa_dd => "JULIA_SDPA_DD_PATH",
                           :sdpa_qd => "JULIA_SDPA_QD_PATH", :sdpa => "JULIA_SDPA_PLAIN_PATH")

"""
    default_use_WSL(variant::Symbol) -> Bool

`true` if the environment variable `JULIA_SDPA_GMP_WSL` (or the `DD`/`QD`/`PLAIN` analogue)
is set to `TRUE`. Then the custom executable from `default_binary_path` is run through WSL.
"""
function default_use_WSL(variant::Symbol)
    var = replace(PATH_ENV_VARS[variant], "_PATH" => "_WSL")
    return uppercase(get(ENV, var, "FALSE")) == "TRUE"
end

"""
    jll_command(variant::Symbol) -> Cmd

The command (with library paths set in its environment) for the executable shipped by the
`_jll` package for `variant`.
"""
function jll_command(variant::Symbol)
    jll = variant == :sdpa_gmp ? SDPA_GMP_jll :
          variant == :sdpa_dd  ? SDPA_DD_jll  :
          variant == :sdpa_qd  ? SDPA_QD_jll  :
          variant == :sdpa     ? SDPA_jll     :
          throw(ArgumentError("unknown variant $variant"))
    if !jll.is_available()
        error("""No prebuilt `$variant` executable is available for this platform
                 ($(Base.BinaryPlatforms.triplet(Base.BinaryPlatforms.HostPlatform()))).
                 Set `$(PATH_ENV_VARS[variant])` to the directory of a `$variant` executable
                 (and `$(replace(PATH_ENV_VARS[variant], "_PATH" => "_WSL"))=TRUE` to run it through WSL),
                 or pass `binary_path` to `SDPAFamily.Optimizer`.""")
    end
    return getproperty(jll, variant)()
end

"""
    WSLize_path(path::String) -> String

This function converts Windows paths for use via WSL.
"""
WSLize_path(path) = replace(path, ":" => "") |> x -> replace(x, "\\" => "/") |> x -> "/mnt/"*x |> lowercase

include("params.jl")
include("objective_function_filter.jl")
include("MOI_wrapper.jl")
include("file_io.jl")
include("presolve.jl")
include("binary_call.jl")

end # module
