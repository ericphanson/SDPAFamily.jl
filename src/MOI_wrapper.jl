# This file modifies code from SDPA.jl (https://github.com/JuliaOpt/SDPA.jl), which is available under an MIT license (see LICENSE).

using MathOptInterface
MOI = MathOptInterface
const MOIU = MOI.Utilities
const AFFEQ{T} = MOI.ConstraintIndex{MOI.ScalarAffineFunction{T}, MOI.EqualTo{T}}

abstract type AbstractBlockMatrix{T} <: AbstractMatrix{T} end
function nblocks end
function block end

struct TemporaryDirectory <: MOI.AbstractOptimizerAttribute end

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

abstract type BlockSolution{T} <: AbstractBlockMatrix{T} end
struct PrimalSolution{T} <: BlockSolution{T}
    blocks::Vector{Matrix{T}}
end
# getptr(X::PrimalSolution, blk) = getResultYMat(X.problem, blk)
struct VarDualSolution{T} <: BlockSolution{T}
    blocks::Vector{Matrix{T}}
end
# getptr(X::VarDualSolution, blk) = getResultXMat(X.problem, blk)
nblocks(X::BlockSolution) = length(X.blocks)
function block(X::BlockSolution, blk::Integer)
    return X.blocks[blk]
end
# Needed by MPB_wrapper
function Base.getindex(A::BlockSolution, i::Integer)
    block(A, i)
end

mutable struct Optimizer{T} <: MOI.AbstractOptimizer
    objconstant::T
    objsign::Int
    blockdims::Vector{Int}
    varmap::Vector{Tuple{Int, Int, Int}} # Variable Index vi -> blk, i, j
    b::Vector{T}
    solve_time::Float64
    verbosity::Verbosity
    options::Dict{Symbol, Any}
    y::Vector{T}
    X::PrimalSolution{T}
    Z::VarDualSolution{T}
    primalobj::T
    dualobj::T
    phasevalue::Symbol
    tempdir::String
    elemdata::Vector{Any}
	presolve::Bool
	binary_path::String
    params::Union{Params, ParamsSetting, String}
    no_solve::Bool
    use_WSL::Bool
	variant::Symbol
    function Optimizer{T}(; variant = :sdpa_gmp, presolve::Bool = false, 
            silent::Bool = false,
            verbose::Verbosity = silent ? SILENT : WARN,
            binary_path = BB_PATHS[variant],
            use_WSL = HAS_WSL[variant],
            params::Union{Params, ParamsSetting, String, NamedTuple} = Params{variant, T}(),
            TemporaryDirectory::String = mktempdir(),
            ) where T

        if params isa NamedTuple
            P = Params{variant, T}(; params...)
        else
            P = params
        end

		optimizer = new(
            zero(T), 1, Int[], Tuple{Int, Int, Int}[], T[], NaN, verbose,
            Dict{Symbol, Any}(), T[], PrimalSolution{T}(Matrix{T}[]),
            VarDualSolution{T}(Matrix{T}[]), zero(T), zero(T), :not_called,
            TemporaryDirectory, [], presolve, binary_path, P, false, use_WSL,
            variant)

        if silent && verbose != SILENT
            throw(ArgumentError("Cannot set both `silent=true` and `verbose != SILENT`."))
        end

		if T != BigFloat && optimizer.verbosity == VERBOSE
			@warn "Not using BigFloat entries may cause underflow errors."
        end

		return optimizer
    end

    Optimizer(; kwargs...) = Optimizer{BigFloat}(; kwargs...)
end

varmap(optimizer::Optimizer, vi::MOI.VariableIndex) = optimizer.varmap[vi.value]

# This code is all broken. TODO: fix.
# function MOI.supports(optimizer::Optimizer, param::MOI.RawParameter)
# 	return param.name in keys(SET_PARAM)
# end
# function MOI.set(optimizer::Optimizer, param::MOI.RawParameter, value)
# 	if !MOI.supports(optimizer, param)
# 		throw(MOI.UnsupportedAttribute(param))
# 	end
# 	optimizer.options[param.name] = value
# end
# function MOI.get(optimizer::Optimizer, param::MOI.RawParameter)
# 	# TODO: This gives a poor error message if the name of the parameter is invalid.
# 	return optimizer.options[param.name]
# end

MOI.supports(::Optimizer, ::MOI.Silent) = true
function MOI.set(optimizer::Optimizer, ::MOI.Silent, value::Bool)
    if value
        optimizer.verbosity = SILENT
    else
        optimizer.verbosity = WARN
    end
end

function MOI.set(optimizer::Optimizer, ::TemporaryDirectory, path::String)
    optimizer.tempdir = path
end

MOI.get(optimizer::Optimizer, ::TemporaryDirectory) = optimizer.tempdir

MOI.get(optimizer::Optimizer, ::MOI.Silent) = optimizer.verbosity == SILENT

MOI.get(::Optimizer, ::MOI.SolverName) = "SDPAFamily"

# See https://www.researchgate.net/publication/247456489_SDPA_SemiDefinite_Programming_Algorithm_User's_Manual_-_Version_600
# "SDPA (SemiDefinite Programming Algorithm) User's Manual — Version 6.00" Section 6.2
const RAW_STATUS = Dict(
    :noINFO => "The iteration has exceeded the maxIteration and stopped with no information on the primal feasibility and the dual feasibility.",
    :pdOPT => "The normal termination yielding both primal and dual approximate optimal solutions.",
    :pFEAS => "The primal problem got feasible but the iteration has exceeded the maxIteration and stopped.",
    :dFEAS => "The dual problem got feasible but the iteration has exceeded the maxIteration and stopped.",
    :pdFEAS => "Both primal problem and the dual problem got feasible, but the iterationhas exceeded the maxIteration and stopped.",
    :pdINF => "At least one of the primal problem and the dual problem is expected to be infeasible.",
    :pFEAS_dINF => "The primal problem has become feasible but the dual problem is expected to be infeasible.",
    :pINF_dFEAS => "The dual problem has become feasible but the primal problem is expected to be infeasible.",
    :pUNBD => "The primal problem is expected to be unbounded.",
    :dUNBD => "The dual problem is expected to be unbounded.")

function MOI.get(optimizer::Optimizer, ::MOI.RawStatusString)
	return RAW_STATUS[optimizer.phasevalue]
end
function MOI.get(optimizer::Optimizer, ::MOI.SolveTime)
	return optimizer.solve_time
end

function MOI.is_empty(optimizer::Optimizer)
    return iszero(optimizer.objconstant) &&
        optimizer.objsign == 1 &&
        isempty(optimizer.blockdims) &&
        isempty(optimizer.varmap) &&
        isempty(optimizer.b) &&
        optimizer.elemdata == []
end
function MOI.empty!(optimizer::Optimizer{T}) where T
    optimizer.objconstant = zero(T)
    optimizer.objsign = 1
    empty!(optimizer.blockdims)
    empty!(optimizer.varmap)
    empty!(optimizer.b)
    optimizer.X = PrimalSolution{T}(Matrix{T}[])
    optimizer.Z = VarDualSolution{T}(Matrix{T}[])
    optimizer.y = T[]
    optimizer.phasevalue = :not_called
    clean_tempdir(optimizer.tempdir)
    optimizer.elemdata = []
    optimizer.primalobj = zero(T)
    optimizer.dualobj = zero(T)
end

function clean_tempdir(tempdir)
    files = joinpath.(Ref(tempdir), ["input.dat-s", "output.dat", "params.sdpa"])
    for f in files
        isfile(f) && rm(f)
    end
end

function MOI.supports(
    optimizer::Optimizer,
    ::Union{MOI.ObjectiveSense,
            MOI.ObjectiveFunction{<:Union{MOI.SingleVariable,
                                          MOI.ScalarAffineFunction{T}}}}) where T
    return true
end

MOI.supports_add_constrained_variables(::Optimizer, ::Type{MOI.Reals}) = false
const SupportedSets = Union{MOI.Nonnegatives, MOI.PositiveSemidefiniteConeTriangle}
MOI.supports_add_constrained_variables(::Optimizer, ::Type{<:SupportedSets}) = true
function MOI.supports_constraint(
    ::Optimizer, ::Type{MOI.VectorOfVariables},
    ::Type{<:SupportedSets})
    return true
end
function MOI.supports_constraint(
    ::Optimizer, ::Type{MOI.ScalarAffineFunction{T}},
    ::Type{MOI.EqualTo{T}}) where T
    return true
end

function MOI.copy_to(dest::Optimizer, src::MOI.ModelLike; kws...)
    return MOIU.automatic_copy_to(dest, src; kws...)
end
MOIU.supports_allocate_load(::Optimizer, copy_names::Bool) = !copy_names

function MOIU.allocate(optimizer::Optimizer, ::MOI.ObjectiveSense, sense::MOI.OptimizationSense)
    # To be sure that it is done before load(optimizer, ::ObjectiveFunction, ...), we do it in allocate
    optimizer.objsign = sense == MOI.MIN_SENSE ? -1 : 1
end
function MOIU.allocate(::Optimizer, ::MOI.ObjectiveFunction, ::Union{MOI.SingleVariable, MOI.ScalarAffineFunction}) end

function MOIU.load(::Optimizer, ::MOI.ObjectiveSense, ::MOI.OptimizationSense) end
# Loads objective coefficient α * vi
function load_objective_term!(optimizer::Optimizer{T}, α, vi::MOI.VariableIndex) where {T}
    blk, i, j = varmap(optimizer, vi)
    coef = optimizer.objsign * α
    if i != j
        coef /= 2
    end
    # in SDP format, it is max and in MPB Conic format it is min
    inputElement(optimizer, 0, blk, i, j, convert(T, coef))
end
function MOIU.load(optimizer::Optimizer, ::MOI.ObjectiveFunction, f::MOI.ScalarAffineFunction)
    obj = MOIU.canonical(f)
    optimizer.objconstant = f.constant
    for t in obj.terms
        if !iszero(t.coefficient)
            load_objective_term!(optimizer, t.coefficient, t.variable_index)
        end
    end
end
function MOIU.load(optimizer::Optimizer{T}, ::MOI.ObjectiveFunction, f::MOI.SingleVariable) where T
    load_objective_term!(optimizer, one(T), f.variable)
end

function new_block(optimizer::Optimizer, set::MOI.Nonnegatives)
    push!(optimizer.blockdims, -MOI.dimension(set))
    blk = length(optimizer.blockdims)
    for i in 1:MOI.dimension(set)
        push!(optimizer.varmap, (blk, i, i))
    end
end

function new_block(optimizer::Optimizer, set::MOI.PositiveSemidefiniteConeTriangle)
    push!(optimizer.blockdims, set.side_dimension)
    blk = length(optimizer.blockdims)
    for i in 1:set.side_dimension
        for j in 1:i
            push!(optimizer.varmap, (blk, i, j))
        end
    end
end

function MOIU.allocate_constrained_variables(optimizer::Optimizer,
                                             set::SupportedSets)
    offset = length(optimizer.varmap)
    new_block(optimizer, set)
    ci = MOI.ConstraintIndex{MOI.VectorOfVariables, typeof(set)}(offset + 1)
    return [MOI.VariableIndex(i) for i in offset .+ (1:MOI.dimension(set))], ci
end

function MOIU.load_constrained_variables(
    optimizer::Optimizer, vis::Vector{MOI.VariableIndex},
    ci::MOI.ConstraintIndex{MOI.VectorOfVariables},
    set::SupportedSets)
end

function MOIU.allocate_variables(model::Optimizer, nvars)
end

function MOIU.load_variables(optimizer::Optimizer{T}, nvars) where T
    @assert nvars == length(optimizer.varmap)
    dummy = isempty(optimizer.b)
    if dummy
        optimizer.b = [one(T)]
        optimizer.blockdims = [optimizer.blockdims; -1]
    end
    if dummy
        inputElement(optimizer, 1, length(optimizer.blockdims), 1, 1, one(T))
    end
end

function MOIU.allocate_constraint(optimizer::Optimizer{T},
                                  func::MOI.ScalarAffineFunction,
                                  set::MOI.EqualTo) where T
    push!(optimizer.b, MOI.constant(set))
    return AFFEQ{T}(length(optimizer.b))
end

function MOIU.load_constraint(m::Optimizer{T}, ci::AFFEQ,
                              f::MOI.ScalarAffineFunction, s::MOI.EqualTo) where T
    if !iszero(MOI.constant(f))
        throw(MOI.ScalarFunctionConstantNotZero{
            T, MOI.ScalarAffineFunction{T}, MOI.EqualTo{T}}(
                MOI.constant(f)))
    end
    f = MOIU.canonical(f) # sum terms with same variables and same outputindex
    for t in f.terms
        if !iszero(t.coefficient)
            blk, i, j = varmap(m, t.variable_index)
            coef = t.coefficient
            if i != j
                coef /= 2
            end
            inputElement(m, ci.value, blk, i, j, convert(T, coef))
        end
    end
end

function MOI.optimize!(m::Optimizer)
	start_time = time()
    # SDPA.initializeUpperTriangle(m.problem, false)
    redundant_F = initializeSolve(m)
    # SDPA.solve(m)
	if m.phasevalue != :pFEAS_dINF
	    inputname = "input.dat-s"
	    outputname = "output.dat"
	    full_input_path = joinpath(m.tempdir, inputname)
	    full_output_path = joinpath(m.tempdir, outputname)
		if !m.no_solve
		    sdpa_gmp_binary_solve!(m, full_input_path, full_output_path, redundant_entries = redundant_F)
		end
	end
    m.solve_time = time() - start_time
end

function MOI.get(m::Optimizer, ::MOI.TerminationStatus)
    status = m.phasevalue
    if status == :not_called
        return MOI.OPTIMIZE_NOT_CALLED
    elseif status == :noINFO
        return MOI.ITERATION_LIMIT
    elseif status == :pFEAS
        return MOI.SLOW_PROGRESS
    elseif status == :dFEAS
        return MOI.SLOW_PROGRESS
    elseif status == :pdFEAS
        return MOI.OPTIMAL
    elseif status == :pdINF
        return MOI.INFEASIBLE_OR_UNBOUNDED
    elseif status == :pFEAS_dINF
        return MOI.INFEASIBLE
    elseif status == :pINF_dFEAS
        return MOI.DUAL_INFEASIBLE
    elseif status == :pdOPT
        return MOI.OPTIMAL
    elseif status == :pUNBD
        return MOI.INFEASIBLE
    elseif status == :dUNBD
        return MOI.DUAL_INFEASIBLE
    end
end

function MOI.get(m::Optimizer, attr::MOI.DualStatus)
    if attr.N > MOI.get(m, MOI.ResultCount())
        return MOI.NO_SOLUTION
    end
    status = m.phasevalue
    if status == :noINFO
        return MOI.UNKNOWN_RESULT_STATUS
    elseif status == :pFEAS
        return MOI.FEASIBLE_POINT
    elseif status == :dFEAS
        return MOI.UNKNOWN_RESULT_STATUS
    elseif status == :pdFEAS
        return MOI.FEASIBLE_POINT
    elseif status == :pdINF
        return MOI.UNKNOWN_RESULT_STATUS
    elseif status == :pFEAS_dINF
        return MOI.INFEASIBILITY_CERTIFICATE
    elseif status == :pINF_dFEAS
        return MOI.INFEASIBLE_POINT
    elseif status == :pdOPT
        return MOI.FEASIBLE_POINT
    elseif status == :pUNBD
        return MOI.INFEASIBILITY_CERTIFICATE
    elseif status == :dUNBD
        return MOI.INFEASIBLE_POINT
    end
end

function MOI.get(m::Optimizer, attr::MOI.PrimalStatus)
    if attr.N > MOI.get(m, MOI.ResultCount())
        return MOI.NO_SOLUTION
    end
    status = m.phasevalue
    if status == :noINFO
        return MOI.UNKNOWN_RESULT_STATUS
    elseif status == :pFEAS
        return MOI.UNKNOWN_RESULT_STATUS
    elseif status == :dFEAS
        return MOI.FEASIBLE_POINT
    elseif status == :pdFEAS
        return MOI.FEASIBLE_POINT
    elseif status == :pdINF
        return MOI.UNKNOWN_RESULT_STATUS
    elseif status == :pFEAS_dINF
        return MOI.INFEASIBLE_POINT
    elseif status == :pINF_dFEAS
        return MOI.INFEASIBILITY_CERTIFICATE
    elseif status == :pdOPT
        return MOI.FEASIBLE_POINT
    elseif status == :pUNBD
        return MOI.INFEASIBLE_POINT
    elseif status == :dUNBD
        return MOI.INFEASIBILITY_CERTIFICATE
    end
end

MOI.get(m::Optimizer, ::MOI.ResultCount) = 1
function MOI.get(m::Optimizer, attr::MOI.ObjectiveValue)
    MOI.check_result_index_bounds(m, attr)
    return m.objsign * m.primalobj + m.objconstant
end

function MOI.get(m::Optimizer, attr::MOI.DualObjectiveValue)
    MOI.check_result_index_bounds(m, attr)
    return m.objsign * m.dualobj + m.objconstant
end
struct PrimalSolutionMatrix <: MOI.AbstractModelAttribute end
MOI.is_set_by_optimize(::PrimalSolutionMatrix) = true
MOI.get(optimizer::Optimizer, ::PrimalSolutionMatrix) = optimizer.X

struct DualSolutionVector <: MOI.AbstractModelAttribute end
MOI.is_set_by_optimize(::DualSolutionVector) = true
function MOI.get(optimizer::Optimizer, ::DualSolutionVector)
    return optimizer.y
end

struct DualSlackMatrix <: MOI.AbstractModelAttribute end
MOI.is_set_by_optimize(::DualSlackMatrix) = true
MOI.get(optimizer::Optimizer, ::DualSlackMatrix) = optimizer.Z

function block(optimizer::Optimizer, ci::MOI.ConstraintIndex{MOI.VectorOfVariables})
    return optimizer.varmap[ci.value][1]
end
function dimension(optimizer::Optimizer, ci::MOI.ConstraintIndex{MOI.VectorOfVariables})
    blockdim = optimizer.blockdims[block(optimizer, ci)]
    if blockdim < 0
        return -blockdim
    else
        return MOI.dimension(MOI.PositiveSemidefiniteConeTriangle(blockdim))
    end
end
function vectorize_block(M, blk::Integer, s::Type{MOI.Nonnegatives})
    return diag(block(M, blk))
end
function vectorize_block(M::AbstractMatrix{T}, blk::Integer, s::Type{MOI.PositiveSemidefiniteConeTriangle}) where T
    B = block(M, blk)
    d = LinearAlgebra.checksquare(B)
    n = MOI.dimension(MOI.PositiveSemidefiniteConeTriangle(d))
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

function MOI.get(optimizer::Optimizer, attr::MOI.VariablePrimal, vi::MOI.VariableIndex)
    MOI.check_result_index_bounds(optimizer, attr)
    blk, i, j = varmap(optimizer, vi)
    return block(MOI.get(optimizer, PrimalSolutionMatrix()), blk)[i, j]
end

function MOI.get(optimizer::Optimizer, attr::MOI.ConstraintPrimal,
                 ci::MOI.ConstraintIndex{MOI.VectorOfVariables, S}) where S<:SupportedSets
    MOI.check_result_index_bounds(optimizer, attr)
    return vectorize_block(MOI.get(optimizer, PrimalSolutionMatrix()), block(optimizer, ci), S)
end

function MOI.get(m::Optimizer, attr::MOI.ConstraintPrimal, ci::AFFEQ)
    MOI.check_result_index_bounds(m, attr)
    return m.b[ci.value]
end

function MOI.get(optimizer::Optimizer, attr::MOI.ConstraintDual,
                 ci::MOI.ConstraintIndex{MOI.VectorOfVariables, S}) where S<:SupportedSets
    MOI.check_result_index_bounds(optimizer, attr)
    return vectorize_block(MOI.get(optimizer, DualSlackMatrix()), block(optimizer, ci), S)
end
function MOI.get(optimizer::Optimizer, attr::MOI.ConstraintDual, ci::AFFEQ)
    MOI.check_result_index_bounds(optimizer, attr)
    return -MOI.get(optimizer, DualSolutionVector())[ci.value]
end
