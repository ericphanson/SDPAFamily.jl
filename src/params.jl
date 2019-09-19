
"""
    struct Params{variant, T <: Number}

An object holding a list of parameters to be used by the solver, parametrized by
the variant and numeric type. `variant` should be one of `:sdpa_gmp, :sdpa_qd,
:sdpa_dd` or `sdpa`. The variant and numeric type are used to choose the
default values for unspecified parameters; they are unused in the case
that every parameter is specified.

!!! note

    It is often simpler to simply pass non-default parameters directly to the
    optimizer as a `NamedTuple`; e.g.
    `SDPAFamily.Optimizer(params = (maxIteration = 600,))`.

# Example

```jldoctest
julia> using SDPAFamily

julia> P = SDPAFamily.Params{:sdpa_gmp, BigFloat}(maxIteration = 600)
SDPAFamily.Params{:sdpa_gmp,BigFloat}(600, 1.0e-30, 10000.0, 2.0, -100000.0, 100000.0, 0.1, 0.3, 0.9, 1.0e-30, 200, "%+.Fe", "%+.Fe", "%+.Fe", "%+.Fe")

julia> SDPAFamily.Optimizer(params = P)
SDPAFamily.Optimizer{BigFloat}

```

# List of parameters

The following is a brief summary of the parameters. See the SDPA manual for more details.

* `maxIteration`: number of iterations allowed * `epsilonStar`: constraint
tolerance * `epsilonDash`: normalized duality gap tolerance * `lambdaStar`:
determines initial point; should have the same order of magnitude as the optimal
solution * `omegaStar`: determines region in which SDPA searches for an optimal
solution; must be at least 1.0. * `lowerBound`, (resp. `upperBound`): bound on
the primal (resp. dual) optimal objective value; serves as a stopping criteria *
`betaStar`: parameter controlling the search direction for feasible points *
`betaBar`: parameter controlling the search direction for infeasible points *
`gammaStar`: reduction factor for the primal and dual step lengths *
`precision`: number of significant bits used for SDPA-GMP; if set to `b` bits,
then `(log(2)/log(10)) * b` is approximately the number of decimal digits of
precision. * `xPrint`, `XPrint`, `YPrint`, `infPrint`: `printf` format
specification used for printing the results to send them from the solver binary
to Julia.

"""
struct Params{variant, T <: Number}
    maxIteration::Int
    epsilonStar::Float64
    lambdaStar::Float64
    omegaStar::Float64
    lowerBound::Float64
    upperBound::Float64
    betaStar::Float64
    betaBar::Float64
    gammaStar::Float64
    epsilonDash::Float64
    precision::Union{Int, Nothing}
    xPrint::Union{String, Nothing}
    XPrint::Union{String, Nothing}
    YPrint::Union{String, Nothing}
    infPrint::Union{String, Nothing}
    function Params{variant, T}(
        maxIteration::Int,
        epsilonStar::Float64,
        lambdaStar::Float64,
        omegaStar::Float64,
        lowerBound::Float64,
        upperBound::Float64,
        betaStar::Float64,
        betaBar::Float64,
        gammaStar::Float64,
        epsilonDash::Float64,
        precision::Union{Int, Nothing} = nothing,
        xPrint::Union{String, Nothing} = nothing,
        XPrint::Union{String, Nothing} = nothing,
        YPrint::Union{String, Nothing} = nothing,
        infPrint::Union{String, Nothing} = nothing) where {variant, T <: Number}

    epsilonStar > 0 || throw(ArgumentError("epsilonStar > 0 violated"))
    lambdaStar > 0 || throw(ArgumentError("lambdaStar > 0 violated"))
    omegaStar > 1 || throw(ArgumentError("omegaStar > 1 violated"))
    1 > betaStar > 0 || throw(ArgumentError("1 > betaStar > 0 violated"))
    1 > betaBar >= 0 || throw(ArgumentError("1 > betaBar >= 0 violated"))
    betaBar >= betaStar || throw(ArgumentError("betaBar >= betaStar violated"))
    1 > gammaStar > 0 || throw(ArgumentError("1 > gammaStar > 0 violated"))
    epsilonDash > 0 || throw(ArgumentError("epsilonDash > 0 violated"))

    if variant == :sdpa_gmp && (precision === nothing || precision < 0)
        throw(ArgumentError("precision > 0 violated for sdpa_gmp"))
    end

    if variant != :sdpa_gmp && precision !== nothing
        throw(ArgumentError("precision should be `nothing` for variants other than `sdpa_gmp`"))
    end

    if any(x -> x === nothing, (xPrint, XPrint, YPrint, infPrint)) != all(x -> x === nothing, (xPrint, XPrint, YPrint, infPrint))
        throw(ArgumentError("Either of xPrint, XPrint, YPrint, infPrint must be specified, or none of them should be."))
    end

        new{variant, T}(maxIteration, epsilonStar, lambdaStar,
        omegaStar, lowerBound, upperBound, betaStar, betaBar,
        gammaStar, epsilonDash, precision, xPrint, XPrint, YPrint, infPrint)
    end
end

function Params{:sdpa_gmp, T}(;
    maxIteration::Int = 400,
    epsilonStar::Float64 = 1e-30,
    lambdaStar::Float64 = 1e4,
    omegaStar::Float64 = 2.0,
    lowerBound::Float64 = -1e5,
    upperBound::Float64 = 1e5,
    betaStar::Float64 = 0.1,
    betaBar::Float64 = 0.3,
    gammaStar::Float64 = 0.9,
    epsilonDash::Float64 = 1e-30,
    precision::Int = 200,
    xPrint::Union{String, Nothing} = "%+.Fe",
    XPrint::Union{String, Nothing} = "%+.Fe",
    YPrint::Union{String, Nothing} = "%+.Fe",
    infPrint::Union{String, Nothing} = "%+.Fe") where {T <: Number}
    Params{:sdpa_gmp, T}(maxIteration, epsilonStar, lambdaStar,
    omegaStar, lowerBound, upperBound, betaStar, betaBar,
    gammaStar, epsilonDash, precision, xPrint, XPrint, YPrint, infPrint)
end

function Params{:sdpa_gmp, Float64}(;
    maxIteration::Int = 200,
    epsilonStar::Float64 = 1e-7,
    lambdaStar::Float64 = 1e3,
    omegaStar::Float64 = 2.0,
    lowerBound::Float64 = -1e5,
    upperBound::Float64 = 1e5,
    betaStar::Float64 = 0.1,
    betaBar::Float64 = 0.3,
    gammaStar::Float64 = 0.9,
    epsilonDash::Float64 = 1e-7,
    precision::Int = 80,
    xPrint::Union{String, Nothing} = "%+.Fe",
    XPrint::Union{String, Nothing} = "%+.Fe",
    YPrint::Union{String, Nothing} = "%+.Fe",
    infPrint::Union{String, Nothing} = "%+.Fe")

    Params{:sdpa_gmp, Float64}(maxIteration, epsilonStar, lambdaStar,
    omegaStar, lowerBound, upperBound, betaStar, betaBar,
    gammaStar, epsilonDash, precision, xPrint, XPrint, YPrint, infPrint)
end

function Params{:sdpa_qd, T}(;
    maxIteration::Int = 100,
    epsilonStar::Float64 = 1e-20,
    lambdaStar::Float64 = 1e3,
    omegaStar::Float64 = 2.0,
    lowerBound::Float64 = -1e5,
    upperBound::Float64 = 1e5,
    betaStar::Float64 = 0.1,
    betaBar::Float64 = 0.3,
    gammaStar::Float64 = 0.9,
    epsilonDash::Float64 = 1e-20) where {T <: Number}

    Params{:sdpa_qd, T}(maxIteration, epsilonStar, lambdaStar,
    omegaStar, lowerBound, upperBound, betaStar, betaBar,
    gammaStar, epsilonDash)
end

function Params{:sdpa_dd, T}(;
    maxIteration::Int = 100,
    epsilonStar::Float64 = 1e-15,
    lambdaStar::Float64 = 1e3,
    omegaStar::Float64 = 2.0,
    lowerBound::Float64 = -1e5,
    upperBound::Float64 = 1e5,
    betaStar::Float64 = 0.1,
    betaBar::Float64 = 0.3,
    gammaStar::Float64 = 0.9,
    epsilonDash::Float64 = 1e-15) where {T <: Number}

    Params{:sdpa_qd, T}(maxIteration, epsilonStar, lambdaStar,
    omegaStar, lowerBound, upperBound, betaStar, betaBar,
    gammaStar, epsilonDash)
end

function Params{:sdpa, T}(;
    maxIteration::Int = 100,
    epsilonStar::Float64 = 1e-7,
    lambdaStar::Float64 = 1e2,
    omegaStar::Float64 = 2.0,
    lowerBound::Float64 = -1e5,
    upperBound::Float64 = 1e5,
    betaStar::Float64 = 0.01,
    betaBar::Float64 = 0.02,
    gammaStar::Float64 = 0.95,
    epsilonDash::Float64 = 1e-7,
    xPrint::Union{String, Nothing} = "%+8.10e",
    XPrint::Union{String, Nothing} = "%+8.10e",
    YPrint::Union{String, Nothing} = "%+8.10e",
    infPrint::Union{String, Nothing} = "%+10.16e") where {T <: Number}

    Params{:sdpa, T}(maxIteration, epsilonStar, lambdaStar,
    omegaStar, lowerBound, upperBound, betaStar, betaBar,
    gammaStar, epsilonDash, nothing, xPrint, XPrint, YPrint, infPrint)
end

function write_params(P::Params{variant}, path::String) where {variant}

    params_string = """$(P.maxIteration)	unsigned int maxIteration;
    $(P.epsilonStar)	double 0.0 < epsilonStar;
    $(P.lambdaStar)   double 0.0 < lambdaStar;
    $(P.omegaStar)   	double 1.0 < omegaStar;
    $(P.lowerBound)  double lowerBound;
    $(P.upperBound)   double upperBound;
    $(P.betaStar)     double 0.0 <= betaStar <  1.0;
    $(P.betaBar)     double 0.0 <= betaBar  <  1.0, betaStar <= betaBar;
    $(P.gammaStar)     double 0.0 < gammaStar  <  1.0;
    $(P.epsilonDash)	double 0.0 < epsilonDash;
    """

    if P.precision !== nothing
        params_string = params_string*"""
        $(P.precision)	precision
        """
    end

    if P.xPrint !== nothing
        params_string = params_string*"""
        $(P.xPrint)     char*  xPrint   (default %+8.3e,   NOPRINT skips printout)
        $(P.XPrint)     char*  XPrint   (default %+8.3e,   NOPRINT skips printout)
        $(P.YPrint)     char*  YPrint   (default %+8.3e,   NOPRINT skips printout)
        $(P.infPrint)     char*  infPrint (default %+10.16e, NOPRINT skips printout)
        """
    end

    open(path, "w") do io
        print(io, params_string)
    end
end
