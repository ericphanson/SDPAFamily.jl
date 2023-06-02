# This file modifies code from SDPA.jl (https://github.com/JuliaOpt/SDPA.jl), which is available under an MIT license (see LICENSE).

import MathOptInterface as MOI
const AFF{T} = MOI.VectorAffineFunction{T}

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
    nn_block_idx::Dict{Int64,Int}
    sd_block_idx::Dict{Int64,Int}
    solve_time::Float64
    verbosity::Verbosity
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
            TemporaryDirectory::String = mktempdir(@get_scratch!("solves")),
            ) where T

        if params isa NamedTuple
            P = Params{variant, T}(; params...)
        else
            P = params
        end

		optimizer = new(
            zero(T), 1, Int[], Dict{Int64,Int}(), Dict{Int64,Int}(), NaN, verbose,
            T[], PrimalSolution{T}(Matrix{T}[]),
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

function MOI.supports(::Optimizer, param::MOI.RawOptimizerAttribute)
	return hasfield(Params, Symbol(param.name))
end
function MOI.set(optimizer::Optimizer, param::MOI.RawOptimizerAttribute, value)
	if !MOI.supports(optimizer, param)
		throw(MOI.UnsupportedAttribute(param))
	end
    setfield!(optimizer.params, Symbol(param.name), value)
    return
end
function MOI.get(optimizer::Optimizer, param::MOI.RawOptimizerAttribute)
	if !MOI.supports(optimizer, param)
		throw(MOI.UnsupportedAttribute(param))
	end
    getfield!(optimizer.params, Symbol(param.name))
    return
end

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
function MOI.get(optimizer::Optimizer, ::MOI.SolveTimeSec)
	return optimizer.solve_time
end

function MOI.is_empty(optimizer::Optimizer)
    return iszero(optimizer.objconstant) &&
        optimizer.objsign == 1 &&
        isempty(optimizer.blockdims) &&
        isempty(optimizer.sd_block_idx) &&
        isempty(optimizer.nn_block_idx) &&
        optimizer.elemdata == []
end
function MOI.empty!(optimizer::Optimizer{T}) where T
    optimizer.objconstant = zero(T)
    optimizer.objsign = 1
    empty!(optimizer.blockdims)
    empty!(optimizer.sd_block_idx)
    empty!(optimizer.nn_block_idx)
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
    ::Optimizer{T},
    ::Union{MOI.ObjectiveSense,
            MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}}) where T
    return true
end

const SupportedSets = Union{MOI.Nonnegatives, MOI.PositiveSemidefiniteConeTriangle}
function MOI.supports_constraint(
    ::Optimizer{T}, ::Type{AFF{T}},
    ::Type{<:SupportedSets}) where {T}
    return true
end

const NONZERO_OBJECTIVE_CONSTANT_ERROR = 
    "Nonzero constant in objective function not supported. Note " *
    "that the constant may be added by the substitution of a " *
    "bridged variable."

function MOI.optimize!(dest::Optimizer{T}, src::MOI.ModelLike) where {T}
    MOI.empty!(dest)
    dest.objsign = MOI.get(src, MOI.ObjectiveSense()) == MOI.MAX_SENSE ? -1 : 1
    inputname = "input.dat-s"
    outputname = "output.dat"
    full_input_path = joinpath(dest.tempdir, inputname)
    full_output_path = joinpath(dest.tempdir, outputname)
    sdpa = MOI.FileFormats.SDPA.Model(number_type = T)
    # SDPA will error if there is a constant so we store it in `objconstant`
    # and filter it out in the copy
    dest.objconstant = MOI.constant(MOI.get(src, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}()))
    index_map = MOI.copy_to(sdpa, ObjectiveFunctionFilter(src))
    # `MOI.FileFormats.SDPA` writes linear block first and then PSD blocks
    NNG = MOI.Nonnegatives
    for ci in MOI.get(sdpa, MOI.ListOfConstraintIndices{AFF{T},NNG}())
        set = MOI.get(sdpa, MOI.ConstraintSet(), ci)
        push!(dest.blockdims, -set.dimension)
        dest.nn_block_idx[ci.value] = length(dest.blockdims)
    end
    PSD = MOI.PositiveSemidefiniteConeTriangle
    for ci in MOI.get(sdpa, MOI.ListOfConstraintIndices{AFF{T},PSD}())
        set = MOI.get(sdpa, MOI.ConstraintSet(), ci)
        push!(dest.blockdims, set.side_dimension)
        dest.sd_block_idx[ci.value] = length(dest.blockdims)
    end
    open(full_input_path, "w") do io
        write(io, sdpa)
    end
	start_time = time()
    if !dest.no_solve
        sdpa_gmp_binary_solve!(dest, full_input_path, full_output_path)
    end
    dest.solve_time = time() - start_time
    return index_map, false
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
        return MOI.DUAL_INFEASIBLE
    elseif status == :pINF_dFEAS
        return MOI.INFEASIBLE
    elseif status == :pdOPT
        return MOI.OPTIMAL
    elseif status == :pUNBD
        return MOI.DUAL_INFEASIBLE
    elseif status == :dUNBD
        return MOI.INFEASIBLE
    end
end

function MOI.get(m::Optimizer, attr::MOI.PrimalStatus)
    if attr.result_index > MOI.get(m, MOI.ResultCount())
        return MOI.NO_SOLUTION
    end
    status = m.phasevalue
    if status == :not_called
        return MOI.NO_SOLUTION
    elseif status == :noINFO
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
    else
        @assert status == :dUNBD
        return MOI.INFEASIBLE_POINT
    end
end

function MOI.get(m::Optimizer, attr::MOI.DualStatus)
    if attr.result_index > MOI.get(m, MOI.ResultCount())
        return MOI.NO_SOLUTION
    end
    status = m.phasevalue
    if status == :not_called
        return MOI.NO_SOLUTION
    elseif status == :noINFO
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
    else
        @assert status == :dUNBD
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

function block(optimizer::Optimizer{T}, ci::MOI.ConstraintIndex{AFF{T},MOI.Nonnegatives}) where {T}
    return optimizer.nn_block_idx[ci.value]
end
function block(optimizer::Optimizer{T}, ci::MOI.ConstraintIndex{AFF{T},MOI.PositiveSemidefiniteConeTriangle}) where {T}
    return optimizer.sd_block_idx[ci.value]
end
function dimension(optimizer::Optimizer, ci::MOI.ConstraintIndex)
    blockdim = optimizer.blockdims[block(optimizer, ci)]
    if blockdim < 0
        return -blockdim
    else
        return MOI.dimension(MOI.PositiveSemidefiniteConeTriangle(blockdim))
    end
end
function vectorize_block(M, blk::Integer, ::Type{MOI.Nonnegatives})
    return diag(block(M, blk))
end
function vectorize_block(M::AbstractMatrix{T}, blk::Integer, ::Type{MOI.PositiveSemidefiniteConeTriangle}) where T
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
    return MOI.get(optimizer, DualSolutionVector())[vi.value]
end

function MOI.get(optimizer::Optimizer, attr::MOI.ConstraintPrimal,
                 ci::MOI.ConstraintIndex{AFF{T}, S}) where {T,S<:SupportedSets}
    MOI.check_result_index_bounds(optimizer, attr)
    return vectorize_block(MOI.get(optimizer, DualSlackMatrix()), block(optimizer, ci), S)
end

function MOI.get(optimizer::Optimizer{T}, attr::MOI.ConstraintDual,
                 ci::MOI.ConstraintIndex{AFF{T}, S}) where {T,S<:SupportedSets}
    MOI.check_result_index_bounds(optimizer, attr)
    return vectorize_block(MOI.get(optimizer, PrimalSolutionMatrix()), block(optimizer, ci), S)
end
