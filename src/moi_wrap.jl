# using MathOptInterface
# MOI = MathOptInterface
# using SemidefiniteOptInterface
# SDOI = SemidefiniteOptInterface

mutable struct SDPAGMPOptimizer{T} <: SDOI.AbstractSDOptimizer
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
    X::SDOI.BlockMatrix{T}
    Z::SDOI.BlockMatrix{T}
    y::Vector{T}
end

SDPAGMPOptimizer{T}() where T = SDPAGMPOptimizer{T}(0,
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
                                                  SDOI.BlockMatrix{T}(Matrix{T}[]),
                                                  SDOI.BlockMatrix{T}(Matrix{T}[]),
                                                  T[])
SDPAGMPoptimizer(T::Type) = SDOI.SDOIOptimizer(SDPAGMPOptimizer{T}(), T)

MOI.get(::SDPAGMPOptimizer, ::MOI.SolverName) = "SDPA_GMP"

function MOI.empty!(mock::SDPAGMPOptimizer{T}) where T
    mock.nconstrs = 0
    mock.blkdims = Int[]
    mock.constraint_constants = T[]
    mock.objective_coefficients = Tuple{T, Int, Int, Int}[]
    mock.constraint_coefficients = Vector{Tuple{T, Int, Int, Int}}[]
    mock.hasprimal = false
    mock.hasdual = false
    mock.terminationstatus = MOI.OPTIMIZE_NOT_CALLED
    mock.primalstatus = MOI.NO_SOLUTION
    mock.dualstatus = MOI.NO_SOLUTION
    mock.X = SDOI.BlockMatrix{T}(Matrix{T}[])
    mock.Z = SDOI.BlockMatrix{T}(Matrix{T}[])
    mock.y = T[]
end
SDOI.coefficienttype(::SDPAGMPOptimizer{T}) where T = T

SDOI.getnumberofconstraints(optimizer::SDPAGMPOptimizer) = optimizer.nconstrs
SDOI.getnumberofblocks(optimizer::SDPAGMPOptimizer) = length(optimizer.blkdims)
SDOI.getblockdimension(optimizer::SDPAGMPOptimizer, blk) = optimizer.blkdims[blk]
function SDOI.init!(optimizer::SDPAGMPOptimizer{T}, blkdims::Vector{Int},
               nconstrs::Integer) where T
    optimizer.nconstrs = nconstrs
    optimizer.blkdims = blkdims
    optimizer.constraint_constants = zeros(T, nconstrs)
    optimizer.objective_coefficients = Tuple{T, Int, Int, Int}[]
    optimizer.constraint_coefficients = map(i -> Tuple{T, Int, Int, Int}[], 1:nconstrs)
end

SDOI.getconstraintconstant(optimizer::SDPAGMPOptimizer, c) = optimizer.constraint_constants[c]
function SDOI.setconstraintconstant!(optimizer::SDPAGMPOptimizer, val, c::Integer)
    optimizer.constraint_constants[c] = val
end

SDOI.getobjectivecoefficients(optimizer::SDPAGMPOptimizer) = optimizer.objective_coefficients
function SDOI.setobjectivecoefficient!(optimizer::SDPAGMPOptimizer, val, blk::Integer, i::Integer, j::Integer)
    push!(optimizer.objective_coefficients, (val, blk, i, j))
end

SDOI.getconstraintcoefficients(optimizer::SDPAGMPOptimizer, c) = optimizer.constraint_coefficients[c]
function SDOI.setconstraintcoefficient!(optimizer::SDPAGMPOptimizer, val, c::Integer,
                                   blk::Integer, i::Integer, j::Integer)
    push!(optimizer.constraint_coefficients[c], (val, blk, i, j))
end

MOI.get(mock::SDPAGMPOptimizer, ::MOI.TerminationStatus) = mock.terminationstatus
function MOI.set(mock::SDPAGMPOptimizer, ::MOI.TerminationStatus,
                 value::MOI.TerminationStatusCode)
    mock.terminationstatus = value
end
MOI.get(mock::SDPAGMPOptimizer, ::MOI.PrimalStatus) = mock.primalstatus
MOI.set(mock::SDPAGMPOptimizer, ::MOI.PrimalStatus, value::MOI.ResultStatusCode) = (mock.primalstatus = value)
MOI.get(mock::SDPAGMPOptimizer, ::MOI.DualStatus) = mock.dualstatus
MOI.set(mock::SDPAGMPOptimizer, ::MOI.DualStatus, value::MOI.ResultStatusCode) = (mock.dualstatus = value)

SDOI.getX(mock::SDPAGMPOptimizer) = mock.X
SDOI.getZ(mock::SDPAGMPOptimizer) = mock.Z
SDOI.gety(mock::SDPAGMPOptimizer) = mock.y

# function MOI.get(m::SDPAGMPOptimizer{T}, ::MOI.DualObjectiveValue) where T
function SDOI.getdualobjectivevalue(mock::SDPAGMPOptimizer{T}) where T
    v = zero(T)
    for (α, blk, i, j) in mock.objective_coefficients
        v += -1*α * block(mock.X, blk)[i, j]
        if i != j
            v += -1*α * block(mock.X, blk)[j, i]
        end
    end
    return v
end
# function MOI.get(m::SDPAGMPOptimizer{T}, ::MOI.ObjectiveValue) where T
function SDOI.getprimalobjectivevalue(mock::SDPAGMPOptimizer{T}) where T
    v = zero(T)
    for c in 1:mock.nconstrs
        v += mock.constraint_constants[c] * mock.y[c]
    end
    return v
end

function MOI.optimize!(m::SDPAGMPOptimizer)
    m.hasprimal = true
    m.hasdual = true
    m.optimize!(m)

    temp = mktempdir()
    inputname = "input.dat-s"
    outputname = "output.dat"
    full_input_path = joinpath(temp, inputname)
    full_output_path = joinpath(temp, outputname)
    MOI.write(m, full_input_path)
    sdpa_gmp_binary_solve(m, full_input_path, full_output_path)
end
