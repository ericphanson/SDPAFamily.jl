abstract type AbstractBlockMatrix{T} <: AbstractMatrix{T} end
struct BlockMatrix{T} <: AbstractBlockMatrix{T}
    blocks::Vector{Matrix{T}}
end
nblocks(bm::BlockMatrix) = length(bm.blocks)
block(bm::BlockMatrix, i::Integer) = bm.blocks[i]

function Base.size(bm::AbstractBlockMatrix)
    n = mapreduce(blk -> LinearAlgebra.checksquare(block(bm, blk)),
                         +, 1:nblocks(bm), init=0)
    return (n, n)
end
function Base.getindex(bm::AbstractBlockMatrix, i::Integer, j::Integer)
    (i < 0 || j < 0) && throw(BoundsError(i, j))
    for k in 1:nblocks(bm)
        blk = block(bm, k)
        n = size(blk, 1)
        if i <= n && j <= n
            return blk[i, j]
        elseif i <= n || j <= n
            return 0
        else
            i -= n
            j -= n
        end
    end
    i, j = (i, j) .+ size(bm)
    throw(BoundsError(i, j))
end
Base.getindex(A::AbstractBlockMatrix, I::Tuple) = getindex(A, I...)


abstract type AbstractSDOptimizer <: MOI.AbstractOptimizer end

const SVF = MOI.SingleVariable
const VVF = MOI.VectorOfVariables
const VF  = Union{SVF, VVF}
const SAF{T} = MOI.ScalarAffineFunction{T}
const ASF{T} = Union{SVF, SAF{T}}

const ZS = Union{MOI.EqualTo, MOI.Zeros}
const NS = Union{MOI.GreaterThan, MOI.Nonnegatives}
const PS = Union{MOI.LessThan, MOI.Nonpositives}
const DS = MOI.PositiveSemidefiniteConeTriangle
const SupportedSets = Union{ZS, NS, PS, DS}

const VI = MOI.VariableIndex
const CI{F, S} = MOI.ConstraintIndex{F, S}

mutable struct SOItoMOIBridge{T, SIT <: AbstractSDOptimizer} <: MOI.AbstractOptimizer
    sdoptimizer::SIT
    setconstant::Dict{Int64, T}
    blkconstant::Dict{Int, T}
    objconstant::T
    objsign::Int
    objshift::T
    nconstrs::Int
    nblocks::Int
    blockdims::Vector{Int}
    free::BitSet
    varmap::Vector{Vector{Tuple{Int, Int, Int, T, T}}} # Variable Index vi -> blk, i, j, coef, shift # x = sum coef * block(X, blk)[i, j] + shift
    zeroblock::Dict{CI, Int}
    constrmap::Dict{CI, UnitRange{Int}} # Constraint Index ci -> cs
    double::Vector{CI} # created when there are two cones for same variable
    function SOItoMOIBridge{T}(sdoptimizer::SIT) where {T, SIT}
        new{T, SIT}(sdoptimizer, Dict{Int64, T}(), Dict{Int, T}(),
            zero(T), 1, zero(T), 0, 0,
            Int[],
            BitSet(),
            Vector{Tuple{Int, Int, Int, T}}[],
            Dict{CI, Int}(),
            Dict{CI, UnitRange{Int}}(),
            CI[])
    end
end
varmap(optimizer::SOItoMOIBridge, vi::VI) = optimizer.varmap[vi.value]
function setvarmap!(optimizer::SOItoMOIBridge{T}, vi::VI, v::Tuple{Int, Int, Int, T, T}) where T
    setvarmap!(optimizer, vi, [v])
end
function setvarmap!(optimizer::SOItoMOIBridge{T}, vi::VI, vs::Vector{Tuple{Int, Int, Int, T, T}}) where T
    optimizer.varmap[vi.value] = vs
end

SDOIOptimizer(sdoptimizer::AbstractSDOptimizer, T=Float64) = SOItoMOIBridge{T}(sdoptimizer)

mutable struct SDPAGMPOptimizer{T} <: AbstractSDOptimizer
    nconstrs::Int
    blkdims::Vector{Int}
    constraint_constants::Vector{T}
    objective_coefficients::Vector{Tuple{T, Int, Int, Int}}
    constraint_coefficients::Vector{Vector{Tuple{T, Int, Int, Int}}}
    optimize!::Function # Function used to set dummy primal/dual values and statuses
    hasprimal::Bool
    hasdual::Bool
    terminationstatus::MOI.TerminationStatusCode
    primalstatus::MOI.ResultStatusCode
    dualstatus::MOI.ResultStatusCode
    X::BlockMatrix{T}
    Z::BlockMatrix{T}
    y::Vector{T}
    verbose::Bool
    normal_sdpa::Bool
end

SDPAGMPOptimizer{T}(verbose::Bool, normal_sdpa::Bool) where T = SDPAGMPOptimizer{T}(0,
                                                  Int[],
                                                  T[],
                                                  Tuple{T, Int, Int, Int}[],
                                                  Vector{Tuple{T, Int, Int, Int}}[],
                                                  (::SDPAGMPOptimizer) -> begin end,
                                                  false,
                                                  false,
                                                  MOI.OPTIMIZE_NOT_CALLED,
                                                  MOI.NO_SOLUTION,
                                                  MOI.NO_SOLUTION,
                                                  BlockMatrix{T}(Matrix{T}[]),
                                                  BlockMatrix{T}(Matrix{T}[]),
                                                  T[],
                                                  verbose,
                                                  normal_sdpa)
SDPAGMPoptimizer(T::Type; verbose = false, normal_sdpa = false)= SDOIOptimizer(SDPAGMPOptimizer{T}(verbose, normal_sdpa), T)

MOI.get(::SDPAGMPOptimizer, ::MOI.SolverName) = "SDPA_GMP"

function MOI.empty!(optimizer::SDPAGMPOptimizer{T}) where T
    optimizer.nconstrs = 0
    optimizer.blkdims = Int[]
    optimizer.constraint_constants = T[]
    optimizer.objective_coefficients = Tuple{T, Int, Int, Int}[]
    optimizer.constraint_coefficients = Vector{Tuple{T, Int, Int, Int}}[]
    optimizer.hasprimal = false
    optimizer.hasdual = false
    optimizer.terminationstatus = MOI.OPTIMIZE_NOT_CALLED
    optimizer.primalstatus = MOI.NO_SOLUTION
    optimizer.dualstatus = MOI.NO_SOLUTION
    optimizer.X = BlockMatrix{T}(Matrix{T}[])
    optimizer.Z = BlockMatrix{T}(Matrix{T}[])
    optimizer.y = T[]
end
coefficienttype(::SDPAGMPOptimizer{T}) where T = T

getnumberofconstraints(optimizer::SDPAGMPOptimizer) = optimizer.nconstrs
getnumberofblocks(optimizer::SDPAGMPOptimizer) = length(optimizer.blkdims)
getblockdimension(optimizer::SDPAGMPOptimizer, blk) = optimizer.blkdims[blk]
function init!(optimizer::SDPAGMPOptimizer{T}, blkdims::Vector{Int},
               nconstrs::Integer) where T
    optimizer.nconstrs = nconstrs
    optimizer.blkdims = blkdims
    optimizer.constraint_constants = zeros(T, nconstrs)
    optimizer.objective_coefficients = Tuple{T, Int, Int, Int}[]
    optimizer.constraint_coefficients = map(i -> Tuple{T, Int, Int, Int}[], 1:nconstrs)
end

getconstraintconstant(optimizer::SDPAGMPOptimizer, c) = optimizer.constraint_constants[c]
function setconstraintconstant!(optimizer::SDPAGMPOptimizer, val, c::Integer)
    optimizer.constraint_constants[c] = val
end

getobjectivecoefficients(optimizer::SDPAGMPOptimizer) = optimizer.objective_coefficients
function setobjectivecoefficient!(optimizer::SDPAGMPOptimizer, val, blk::Integer, i::Integer, j::Integer)
    push!(optimizer.objective_coefficients, (val, blk, i, j))
end

getconstraintcoefficients(optimizer::SDPAGMPOptimizer, c) = optimizer.constraint_coefficients[c]
function setconstraintcoefficient!(optimizer::SDPAGMPOptimizer, val, c::Integer,
                                   blk::Integer, i::Integer, j::Integer)
    push!(optimizer.constraint_coefficients[c], (val, blk, i, j))
end

MOI.get(optimizer::SDPAGMPOptimizer, ::MOI.TerminationStatus) = optimizer.terminationstatus
function MOI.set(optimizer::SDPAGMPOptimizer, ::MOI.TerminationStatus,
                 value::MOI.TerminationStatusCode)
    optimizer.terminationstatus = value
end
MOI.get(optimizer::SDPAGMPOptimizer, ::MOI.PrimalStatus) = optimizer.primalstatus
MOI.set(optimizer::SDPAGMPOptimizer, ::MOI.PrimalStatus, value::MOI.ResultStatusCode) = (optimizer.primalstatus = value)
MOI.get(optimizer::SDPAGMPOptimizer, ::MOI.DualStatus) = optimizer.dualstatus
MOI.set(optimizer::SDPAGMPOptimizer, ::MOI.DualStatus, value::MOI.ResultStatusCode) = (optimizer.dualstatus = value)

getX(optimizer::SDPAGMPOptimizer) = optimizer.X
getZ(optimizer::SDPAGMPOptimizer) = optimizer.Z
gety(optimizer::SDPAGMPOptimizer) = optimizer.y

# function MOI.get(m::SDPAGMPOptimizer{T}, ::MOI.DualObjectiveValue) where T
function getdualobjectivevalue(optimizer::SDPAGMPOptimizer{T}) where T
    v = zero(T)
    for (α, blk, i, j) in optimizer.objective_coefficients
        v += -1*α * block(optimizer.X, blk)[i, j]
        if i != j
            v += -1*α * block(optimizer.X, blk)[j, i]
        end
    end
    return v
end
# function MOI.get(m::SDPAGMPOptimizer{T}, ::MOI.ObjectiveValue) where T
function getprimalobjectivevalue(optimizer::SDPAGMPOptimizer{T}) where T
    v = zero(T)
    for c in 1:optimizer.nconstrs
        v += optimizer.constraint_constants[c] * optimizer.y[c]
    end
    return v
end

function MOI.optimize!(optimizer::SDPAGMPOptimizer{T}) where T
    optimizer.hasprimal = true
    optimizer.hasdual = true
    optimizer.optimize!(optimizer)

    temp = mktempdir()
    inputname = "input.dat-s"
    outputname = "output.dat"
    full_input_path = joinpath(temp, inputname)
    full_output_path = joinpath(temp, outputname)
    MOI.write(optimizer, full_input_path)
    sdpa_gmp_binary_solve!(optimizer, full_input_path, full_output_path)
    if optimizer.verbose
        for i in readlines(full_output_path)
            println(stdout, i)
        end
    end
end



"""
code from SDOI
"""


# include("load.jl")

# include("variable.jl")
const VIS = Union{VI, Vector{VI}}

function newblock(m::SOItoMOIBridge, n)
    push!(m.blockdims, n)
    m.nblocks += 1
end

isfree(m, v::VI) = v.value in m.free
function unfree(m, v)
    @assert isfree(m, v)
    delete!(m.free, v.value)
end

function _constraintvariable!(m::SOItoMOIBridge{T}, vs::VIS, s::ZS) where T
    blk = newblock(m, -_length(vs))
    for (i, v) in _enumerate(vs)
        setvarmap!(m, v, (blk, i, i, one(T), _getconstant(m, s)))
        unfree(m, v)
    end
    blk
end
vscaling(::Type{<:NS}, T) = one(T) #1.
vscaling(::Type{<:PS}, T) = -1*one(T)#-1.
_length(vi::VI) = 1
_length(vi::Vector{VI}) = length(vi)
_enumerate(vi::VI) = enumerate((vi,))
_enumerate(vi::Vector{VI}) = enumerate(vi)
function _constraintvariable!(m::SOItoMOIBridge{T}, vs::VIS, s::S) where {S<:Union{NS, PS}, T}
    blk = newblock(m, -_length(vs))
    cst = _getconstant(m, s)
    m.blkconstant[blk] = cst
    for (i, v) in _enumerate(vs)
        setvarmap!(m, v, (blk, i, i, vscaling(S, T), cst))
        unfree(m, v)
    end
    blk
end
function getmatdim(k::Integer)
    # n*(n+1)/2 = k
    # n^2+n-2k = 0
    # (-1 + sqrt(1 + 8k))/2
    n = div(isqrt(1 + 8k) - 1, 2)
    if n * (n+1) != 2*k
        error("sd dim not consistent")
    end
    n
end
function _constraintvariable!(m::SOItoMOIBridge{T}, vs::VIS, ::DS) where T
    d = getmatdim(length(vs))
    k = 0
    blk = newblock(m, d)
    for i in 1:d
        for j in 1:i
            k += 1
            setvarmap!(m, vs[k], (blk, i, j, i == j ? one(T) : one(T)/2, zero(T)))
            unfree(m, vs[k])
        end
    end
    blk
end
_var(f::SVF) = f.variable
_var(f::VVF) = f.variables
function _throw_error_if_unfree(m, vi::MOI.VariableIndex)
    if !isfree(m, vi)
        error("A variable cannot be constrained by multiple ",
              "`MOI.SingleVariable` or `MOI.VectorOfVariables` constraints.")
    end
end
function _throw_error_if_unfree(m, vis::MOI.Vector)
    for vi in vis
        _throw_error_if_unfree(m, vi)
    end
end
function MOIU.allocate_constraint(m::SOItoMOIBridge{T}, f::VF, s::SupportedSets) where T
    vis = _var(f)
    _throw_error_if_unfree(m, vis)
    blk = _constraintvariable!(m, vis, s)
    if isa(s, ZS)
        ci = _allocate_constraint(m, f, s)
        m.zeroblock[ci] = blk
        return ci
    else
        return CI{typeof(f), typeof(s)}(-blk)
    end
end

_getconstant(m::SOItoMOIBridge, s::MOI.AbstractScalarSet) = MOIU.getconstant(s)
_getconstant(m::SOItoMOIBridge{T}, s::MOI.AbstractSet) where T = zero(T)

_var(f::SVF, j) = f.variable
_var(f::VVF, j) = f.variables[j]
function MOIU.load_constraint(m::SOItoMOIBridge, ci::CI, f::VF, s::SupportedSets)
    if ci.value >= 0 # i.e. s is ZS or _var(f) wasn't free at allocate_constraint
        setconstant!(m, ci, s)
        cs = m.constrmap[ci]
        @assert !isempty(cs)
        for k in 1:length(cs)
            vm = varmap(m, _var(f, k))
            # For free variables, the length of vm is 2, clearly not the case here
            @assert length(vm) == 1
            (blk, i, j, coef, shift) = first(vm)
            c = cs[k]
            setconstraintcoefficient!(m.sdoptimizer, coef, c, blk, i, j)
            setconstraintconstant!(m.sdoptimizer,  _getconstant(m, s) - coef * shift, c)
        end
    end
end

function loadfreevariables!(m::SOItoMOIBridge{T}) where T
    for vi in m.free
        blk = newblock(m, -2)
        # x free transformed into x = y - z with y, z >= 0
        setvarmap!(m, VI(vi), [(blk, 1, 1, one(T), zero(T)), (blk, 2, 2, -one(T), zero(T))])
    end
end



# include("constraint.jl")

nconstraints(f::Union{SVF, SAF}, s) = 1
nconstraints(f::VVF, s) = length(f.variables)

function _allocate_constraint(m::SOItoMOIBridge, f, s)
    ci = CI{typeof(f), typeof(s)}(m.nconstrs)
    n = nconstraints(f, s)
    # Fails on Julia v0.6
    #m.constrmap[ci] = m.nconstrs .+ (1:n)
    m.constrmap[ci] = (m.nconstrs + 1):(m.nconstrs + n)
    m.nconstrs += n
    return ci
end
function MOIU.allocate_constraint(m::SOItoMOIBridge, f::SAF, s::SupportedSets)
    _allocate_constraint(m::SOItoMOIBridge, f, s)
end

function loadcoefficients!(m::SOItoMOIBridge, cs::UnitRange,
                           f::MOI.ScalarAffineFunction, s)
    f = MOIU.canonical(f) # sum terms with same variables and same outputindex
    @assert length(cs) == 1
    c = first(cs)
    rhs = MOIU.getconstant(s) - MOI._constant(f)
    for t in f.terms
        if !iszero(t.coefficient)
            for (blk, i, j, coef, shift) in varmap(m, t.variable_index)
                if !iszero(blk)
                    @assert !iszero(coef)
                    setconstraintcoefficient!(m.sdoptimizer, t.coefficient*coef, c, blk, i, j)
                end
                rhs -= t.coefficient * shift
            end
        end
    end
    setconstraintconstant!(m.sdoptimizer, rhs, c)
end

function MOIU.load_constraint(m::SOItoMOIBridge, ci::CI, f::SAF, s::SupportedSets)
    setconstant!(m, ci, s)
    cs = m.constrmap[ci]
    @assert !isempty(cs)
    loadcoefficients!(m, cs, f, s)
end

# load.jl

function MOIU.allocate(optimizer::SOItoMOIBridge, ::MOI.ObjectiveSense, sense::MOI.OptimizationSense)
    # To be sure that it is done before load(optimizer, ::ObjectiveFunction, ...), we do it in allocate
    optimizer.objsign = sense == MOI.MIN_SENSE ? -1 : 1
end
function MOIU.allocate(::SOItoMOIBridge, ::MOI.ObjectiveFunction, ::Union{MOI.SingleVariable, MOI.ScalarAffineFunction}) end

function MOIU.load(::SOItoMOIBridge, ::MOI.ObjectiveSense, ::MOI.OptimizationSense) end
# Loads objective coefficient α * vi
function load_objective_term!(optimizer::SOItoMOIBridge, α, vi::MOI.VariableIndex)
    for (blk, i, j, coef, shift) in varmap(optimizer, vi)
        if !iszero(blk)
            # in SDP format, it is max and in MPB Conic format it is min
            setobjectivecoefficient!(optimizer.sdoptimizer, optimizer.objsign * coef * α, blk, i, j)
        end
        optimizer.objshift += α * shift
    end
end
function MOIU.load(optimizer::SOItoMOIBridge, ::MOI.ObjectiveFunction, f::MOI.ScalarAffineFunction)
    obj = MOIU.canonical(f)
    optimizer.objconstant = f.constant
    for t in obj.terms
        if !iszero(t.coefficient)
            load_objective_term!(optimizer, t.coefficient, t.variable_index)
        end
    end
end
function MOIU.load(optimizer::SOItoMOIBridge{T}, ::MOI.ObjectiveFunction, f::MOI.SingleVariable) where T
    load_objective_term!(optimizer, one(T), f.variable)
end

function MOIU.allocate_variables(optimizer::SOItoMOIBridge{T}, nvars) where T
    optimizer.free = BitSet(1:nvars)
    optimizer.varmap = Vector{Vector{Tuple{Int, Int, Int, T, T}}}(undef, nvars)
    VI.(1:nvars)
end

function MOIU.load_variables(optimizer::SOItoMOIBridge, nvars)
    @assert nvars == length(optimizer.varmap)
    loadfreevariables!(optimizer)
    init!(optimizer.sdoptimizer, optimizer.blockdims, optimizer.nconstrs)
end


function MOI.get(optimizer::SOItoMOIBridge, attr::MOI.SolverName)
    return MOI.get(optimizer.sdoptimizer, attr)
end

function MOI.is_empty(optimizer::SOItoMOIBridge)
    isempty(optimizer.double) &&
    isempty(optimizer.setconstant) &&
    isempty(optimizer.blkconstant) &&
    iszero(optimizer.objconstant) &&
    optimizer.objsign == 1 &&
    iszero(optimizer.objshift) &&
    iszero(optimizer.nconstrs) &&
    iszero(optimizer.nblocks) &&
    isempty(optimizer.blockdims) &&
    isempty(optimizer.free) &&
    isempty(optimizer.varmap) &&
    isempty(optimizer.zeroblock) &&
    isempty(optimizer.constrmap)
end
function MOI.empty!(optimizer::SOItoMOIBridge{T}) where T
    for s in optimizer.double
        MOI.delete(optimizer, s)
    end
    MOI.empty!(optimizer.sdoptimizer)
    optimizer.double = CI[]
    optimizer.setconstant = Dict{Int64, T}()
    optimizer.blkconstant = Dict{Int, T}()
    optimizer.objconstant = zero(T)
    optimizer.objsign = 1
    optimizer.objshift = zero(T)
    optimizer.nconstrs = 0
    optimizer.nblocks = 0
    optimizer.blockdims = Int[]
    optimizer.free = BitSet()
    optimizer.varmap = Vector{Tuple{Int, Int, Int, T}}[]
    optimizer.zeroblock = Dict{CI, Int}()
    optimizer.constrmap = Dict{CI, UnitRange{Int}}()
end

function setconstant!(optimizer::SOItoMOIBridge, ci::CI, s) end
function setconstant!(optimizer::SOItoMOIBridge, ci::CI, s::MOI.AbstractScalarSet)
    optimizer.setconstant[ci.value] = MOIU.getconstant(s)
end
function set_constant(optimizer::SOItoMOIBridge,
                      ci::CI{<:MOI.AbstractScalarFunction,
                             <:MOI.AbstractScalarSet})
    return optimizer.setconstant[ci.value]
end
function set_constant(optimizer::SOItoMOIBridge{T}, ci::CI) where T
    return zeros(T, length(optimizer.constrmap[ci]))
end
function addblkconstant(optimizer::SOItoMOIBridge, ci::CI{<:Any, <:Union{NS, PS}}, x)
    blk = -ci.value
    return x .+ optimizer.blkconstant[blk]
end
addblkconstant(optimizer::SOItoMOIBridge, ci::CI, x) = x

function MOI.supports(
    optimizer::SOItoMOIBridge{T},
    ::Union{MOI.ObjectiveSense,
            MOI.ObjectiveFunction{<:Union{MOI.SingleVariable,
                                          MOI.ScalarAffineFunction{T}}}}) where T
    return true
end

# Zeros and Nonpositives supports could be removed thanks to variable bridges
# * `VectorOfVariables`-in-`Zeros` would return a `VectorAffineFunction` with
#   zero constant and no variable created.
# * `VectorOfVariables`-in-`Nonpositives` would create variables in
#   `Nonnegatives` and return a `VectorAffineFunction` containing `-` the
#    variables.
function MOI.supports_constraint(
    ::SOItoMOIBridge, ::Type{MOI.VectorOfVariables},
    ::Type{<:Union{MOI.Zeros, MOI.Nonnegatives, MOI.Nonpositives,
                   MOI.PositiveSemidefiniteConeTriangle}})
    return true
end
# This support could be remove thanks to variable bridges.
# The VectorizeVariableBridge would redirect to the above case and then the
# resulting function would be shifted by the constant.
function MOI.supports_constraint(
    ::SOItoMOIBridge{T}, ::Type{MOI.SingleVariable},
    ::Type{<:Union{MOI.EqualTo{T}, MOI.GreaterThan{T}, MOI.LessThan{T}}}) where T
    return true
end
function MOI.supports_constraint(
    ::SOItoMOIBridge{T}, ::Type{MOI.ScalarAffineFunction{T}},
    ::Type{MOI.EqualTo{T}}) where T
    return true
end

function MOI.copy_to(dest::SOItoMOIBridge, src::MOI.ModelLike; kws...)
    return MOIU.automatic_copy_to(dest, src; kws...)
end
MOIU.supports_allocate_load(::SOItoMOIBridge, copy_names::Bool) = !copy_names

MOI.optimize!(m::SOItoMOIBridge) = MOI.optimize!(m.sdoptimizer)

# Objective

function MOI.get(m::SOItoMOIBridge, ::MOI.ObjectiveValue)
    m.objshift + m.objsign * getprimalobjectivevalue(m.sdoptimizer) + m.objconstant
end

# Attributes

const SolverStatus = Union{MOI.TerminationStatus, MOI.PrimalStatus, MOI.DualStatus}
MOI.get(m::SOItoMOIBridge, s::SolverStatus) = MOI.get(m.sdoptimizer, s)


MOI.get(m::SOItoMOIBridge, ::MOI.ResultCount) = 1

function _getblock(M, blk::Integer, s::Type{<:Union{NS, ZS}})
    return diag(block(M, blk))
end
function _getblock(M, blk::Integer, s::Type{<:PS})
    return -diag(block(M, blk))
end
# Vectorized length for matrix dimension d
sympackedlen(d::Integer) = (d*(d+1)) >> 1
function _getblock(M::AbstractMatrix{T}, blk::Integer, s::Type{<:DS}) where T
    B = block(M, blk)
    d = LinearAlgebra.checksquare(B)
    n = sympackedlen(d)
    v = Vector{T}(undef, n)
    k = 0
    for j in 1:d
        for i in 1:j
            k += 1
            v[k] = B[i, j]
        end
    end
    @assert k == n
    return v
end
function getblock(M, blk::Integer, s::Type{<:MOI.AbstractScalarSet})
    vd = _getblock(M, blk, s)
    @assert length(vd) == 1
    return vd[1]
end
function getblock(M, blk::Integer, s::Type{<:MOI.AbstractVectorSet})
    return _getblock(M, blk, s)
end

getvarprimal(m::SOItoMOIBridge, blk::Integer, S) = getblock(getX(m.sdoptimizer), blk, S)
function getvardual(m::SOItoMOIBridge, blk::Integer, S)
    z = getZ(m.sdoptimizer)
    b = getblock(z, blk, S)
    return getblock(getZ(m.sdoptimizer), blk, S)
end

function MOI.get(m::SOItoMOIBridge{T}, ::MOI.VariablePrimal, vi::VI) where T
    X = getX(m.sdoptimizer)
    x = zero(T)
    for (blk, i, j, coef, shift) in varmap(m, vi)
        x += shift
        if blk != 0
            x += block(X, blk)[i, j] * sign(coef)
        end
    end
    return x
end
function MOI.get(m::SOItoMOIBridge, vp::MOI.VariablePrimal, vi::Vector{VI})
    return MOI.get.(m, vp, vi)
end

function _getattribute(m::SOItoMOIBridge, ci::CI{<:ASF}, f)
    cs = m.constrmap[ci]
    @assert length(cs) == 1
    return f(m, first(cs))
end
function _getattribute(m::SOItoMOIBridge, ci::CI{<:VVF}, f)
    return f.(m, m.constrmap[ci])
end

function MOI.get(m::SOItoMOIBridge, a::MOI.ConstraintPrimal,
                 ci::CI{F, S}) where {F, S}
    if ci.value >= 0
        return set_constant(m, ci)
    else
        # Variable Function-in-S with S different from Zeros and EqualTo and not a double variable constraint
        blk = -ci.value
        return addblkconstant(m, ci, getvarprimal(m, blk, S))
    end
end

function MOI.get(m::SOItoMOIBridge, ::MOI.ConstraintDual, ci::CI{<:VF, S}) where S<:SupportedSets
    if ci.value < 0
        return getvardual(m, -ci.value, S)
    else
        dual = _getattribute(m, ci, getdual)
        if haskey(m.zeroblock, ci) # ZS
            return dual + getvardual(m, m.zeroblock[ci], S)
        else # var constraint on unfree constraint
            return dual
        end
    end
end

function getdual(m::SOItoMOIBridge{T}, c::Integer) where T
    if c == 0
        return zero(T)
    else
        return -gety(m.sdoptimizer)[c]
    end
end
function MOI.get(m::SOItoMOIBridge, ::MOI.ConstraintDual, ci::CI)
    return _getattribute(m, ci, getdual)
end

# include("sdpa.jl")

# SDPA file format reader-writer

# Redefine the following method to branch on the filename extension of `optimizer` supports more formats
MOI.write(optimizer::SOItoMOIBridge, filename::String) = MOI.write(optimizer.sdoptimizer, filename)
MOI.write(optimizer::AbstractSDOptimizer, filename::String) = writesdpa(optimizer, filename)
MOI.read!(optimizer::SOItoMOIBridge, filename::String) = MOI.read!(optimizer.sdoptimizer, filename)
MOI.read!(optimizer::AbstractSDOptimizer, filename::String) = readsdpa!(optimizer, filename)

function writesdpa(optimizer::AbstractSDOptimizer, filename::String)
    endswith(filename, ".dat-s") || @warn "Filename must end with .dat-s $filename"
    file = open(filename, "w") do io
        nconstrs = getnumberofconstraints(optimizer)
        println(io, nconstrs)
        nblocks = getnumberofblocks(optimizer)
        println(io, nblocks)
        for blk in 1:nblocks
            print(io, getblockdimension(optimizer, blk))
            if blk != nblocks
                print(io, ' ')
            end
        end
        println(io)
        for c in 1:nconstrs
            print(io, getconstraintconstant(optimizer, c))
            if c != nconstrs
                print(io, ' ')
            end
        end
        println(io)
        for (val, blk, i, j) in getobjectivecoefficients(optimizer)
            val = val
            println(io, "0 $blk $i $j $val")
        end
        for c in 1:nconstrs
            for (val, blk, i, j) in getconstraintcoefficients(optimizer, c)
                val = val
                println(io, "$c $blk $i $j $val")
            end
        end
    end
end

nextline(io::IO) = chomp(readline(io))

function readsdpa!(optimizer::AbstractSDOptimizer, filename::String)
    endswith(filename, ".dat-s") || error("Filename '$filename' must end with .dat-s")
    open(filename, "r") do io
        line = nextline(io)
        while line[1] == '"' || line[1] == '*' # Comment
            line = nextline(io)
        end
        nconstrs = parse(Int, line)
        nblocks = parse(Int, nextline(io))
        blkdims = parse.(Int, split(nextline(io)))
        init!(optimizer, blkdims, nconstrs)
        T = coefficienttype(optimizer)
        constraint_constants = parse.(T, split(nextline(io)))
        for c in 1:nconstrs
            setconstraintconstant!(optimizer, constraint_constants[c], c)
        end
        while !eof(io)
            line = nextline(io)
            isempty(line) && break
            s = split(line)
            c = parse(Int, s[1])
            0 ≤ c ≤ nconstrs || error("Invalid constraint index $c in '$filename', it should be an integer between 0 and $nconstrs")
            blk = parse(Int, s[2])
            1 ≤ blk ≤ nblocks || error("Invalid block index $blk in '$filename', it should be an integer between 0 and $nblocks")
            i = parse(Int, s[3])
            j = parse(Int, s[4])
            val = parse(T, s[5])
            if iszero(c)
                setobjectivecoefficient!(optimizer, val, blk, i, j)
            else
                setconstraintcoefficient!(optimizer, val, c, blk, i, j)
            end
        end
    end
end