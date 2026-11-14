using SparseArrays

# # Presolve
#
# `SDPAFamily.Optimizer` receives a model of the form
#
#     minimize    c'x
#     subject to  A_K x + b_K ∈ K        (`Nonnegatives` and `PositiveSemidefiniteConeTriangle`)
#                 A_E x + b_E = 0        (`Zeros`; only supported with `presolve = true`)
#
# with free variables `x`. SDPA solves the pair of problems
#
#     (P)  minimize   c'x          subject to  F_1 x_1 + ... + F_m x_m - F_0 ⪰ 0
#     (D)  maximize   F_0 • Y      subject to  F_i • Y = c_i for i = 1, ..., m,  Y ⪰ 0,
#
# where the matrix `F_i` is column `i` of `A_K`. Two properties of the model cause
# the primal-dual interior point method of SDPA to fail or to return inaccurate results:
#
# 1. Equality constraints. Without presolve, MathOptInterface bridges `A_E x + b_E = 0`
#    into the two inequalities `A_E x + b_E ≥ 0` and `-A_E x - b_E ≥ 0`. Then (P) has no
#    strictly feasible point.
# 2. Linearly dependent columns of `A_K`. Then the equality constraints of (D) are
#    linearly dependent, and the Schur complement matrix that SDPA factorizes in each
#    iteration is singular.
#
# The presolve removes both, in the number type `T` of the optimizer:
#
# 1. Gaussian elimination on the rows of `A_E` selects one pivot variable per linearly
#    independent row and expresses it as an affine function of the remaining ("free")
#    variables. Substituting these expressions into the cone constraints and the
#    objective gives a model without equality constraints. Linearly dependent rows of
#    `A_E` are dropped if they are consistent, and otherwise the model is infeasible
#    and a dual certificate is returned.
# 2. Gaussian elimination on the columns of the substituted `A_K` finds a maximal set
#    of linearly independent columns. The other variables are fixed to zero. If
#    `F_k = Σᵢ aᵢ F_i` and `c_k = Σᵢ aᵢ c_i`, this changes neither the set of values of
#    `F_1 x_1 + ... + F_m x_m` nor the objective values that can be reached, and the
#    equality `F_k • Y = c_k` of (D) is implied by the others. If `c_k ≠ Σᵢ aᵢ c_i`,
#    then (D) is infeasible, and `d = e_k - Σᵢ aᵢ e_i` (with the sign such that
#    `c'd < 0`) is a primal ray that certifies this.
#
# After SDPA has solved the reduced problem, the solution is mapped back:
#
# * The values of the eliminated variables are computed from the free variables.
# * The primal and dual values of the cone constraints are the same in the reduced and
#   the original model, since the reduced constraints have the same affine values and
#   the dual matrix `Y` satisfies the dual constraints of both models.
# * The duals `y_E` of the equality constraints solve `A_E' y_E = c - A_K' y_K`. The
#   right-hand side is in the row space of `A_E` when `y_K` is dual feasible for the
#   reduced model, and the system is solved with the same elimination. Dropped
#   (dependent) rows get the dual value zero.

default_presolve_tol(::Type{T}) where {T} = T <: Union{Rational,Integer} ? zero(T) : sqrt(eps(T))
_drop_tol(::Type{T}) where {T} = T <: Union{Rational,Integer} ? zero(T) : 8 * eps(T)

"""
    Elimination{T}

Result of [`SDPAFamily.eliminate`](@ref).
"""
mutable struct Elimination{T}
    # basis vector `q` has the entry `one(T)` in row `pivots[q]`, and zeros in the pivot
    # rows of the basis vectors before it
    pivots::Vector{Int}
    basis::Vector{Dict{Int,T}}
    basis_c::Vector{T}
    # basis vector `q` is `Σ_j comb[q][j] * A[:, j]` (only with `track = true`)
    comb::Vector{Dict{Int,T}}
    independent::Vector{Int}
    redundant::Vector{Int}
    # the first column that is dependent but inconsistent with `c`, or 0
    inconsistent::Int
    # with `track = true` and `inconsistent > 0`: coefficients `d` with `A d ≈ 0`, `c'd < 0`
    certificate::Dict{Int,T}
end

"""
    eliminate(A::SparseMatrixCSC{T}, c::AbstractVector{T}; tol = default_presolve_tol(T), track = false) -> Elimination{T}

Gaussian elimination with partial pivoting on the columns of `A`, in the number type
`T`. The columns are processed in order. Column `k` is linearly dependent on the
previous columns if, after elimination, its largest remaining entry is at most `tol`
times its largest original entry. For such a column the same combination is applied to
`c`; if the remainder of `c[k]` is larger than `tol` times a bound on the size of the
combined terms, column `k` is recorded as inconsistent and the elimination stops.

With `tol = 0` (the default for `Rational` and `Integer` element types) the
computation is exact.
"""
function eliminate(A::SparseMatrixCSC{T}, c::AbstractVector{T};
                   tol = default_presolve_tol(T), track::Bool = false) where {T}
    n = size(A, 2)
    length(c) == n || throw(DimensionMismatch("`c` must have one entry per column of `A`"))
    cscale = isempty(c) ? zero(T) : maximum(abs, c)
    E = Elimination{T}(Int[], Dict{Int,T}[], T[], Dict{Int,T}[], Int[], Int[], 0,
                       Dict{Int,T}())
    # an upper bound on the 1-norm of the coefficients of each basis vector in terms of
    # the original columns; it scales the tolerance of the consistency check
    combsize = T[]
    row_to_basis = Dict{Int,Int}()
    rows_A = rowvals(A)
    vals_A = nonzeros(A)
    for k in 1:n
        v = Dict{Int,T}()
        for p in nzrange(A, k)
            iszero(vals_A[p]) && continue
            v[rows_A[p]] = get(v, rows_A[p], zero(T)) + vals_A[p]
        end
        scale = maximum(abs, values(v); init = zero(T))
        droptol = _drop_tol(T) * scale
        cv = c[k]
        comb = track ? Dict{Int,T}(k => one(T)) : Dict{Int,T}()
        size_k = one(T)
        todo = BitSet()
        for r in keys(v)
            q = get(row_to_basis, r, 0)
            q > 0 && push!(todo, q)
        end
        # Eliminate with the basis vectors in the order they were added: basis vector
        # `q` is zero in the pivot rows of all basis vectors added before it.
        while !isempty(todo)
            q = first(todo)
            delete!(todo, q)
            prow = E.pivots[q]
            a = get(v, prow, zero(T))
            iszero(a) && continue
            for (r, val) in E.basis[q]
                r == prow && continue
                newval = get(v, r, zero(T)) - a * val
                if abs(newval) <= droptol
                    delete!(v, r)
                else
                    v[r] = newval
                    q2 = get(row_to_basis, r, 0)
                    q2 > q && push!(todo, q2)
                end
            end
            delete!(v, prow)
            cv -= a * E.basis_c[q]
            size_k += abs(a) * combsize[q]
            if track
                for (j, w) in E.comb[q]
                    comb[j] = get(comb, j, zero(T)) - a * w
                end
            end
        end
        resid = maximum(abs, values(v); init = zero(T))
        if resid <= tol * scale
            # dependent column (this includes a column of zeros)
            if abs(cv) > tol * cscale * size_k
                E.inconsistent = k
                if track
                    s = cv > 0 ? -one(T) : one(T)
                    for (j, w) in comb
                        iszero(w) || (E.certificate[j] = s * w)
                    end
                end
                return E
            end
            push!(E.redundant, k)
        else
            # new basis vector; the pivot is the largest entry (smallest row on ties)
            prow, pval = 0, zero(T)
            for (r, val) in v
                if abs(val) > abs(pval) || (abs(val) == abs(pval) && r < prow)
                    prow, pval = r, val
                end
            end
            push!(E.pivots, prow)
            push!(E.basis, Dict{Int,T}(r => val / pval for (r, val) in v))
            E.basis[end][prow] = one(T)
            push!(E.basis_c, cv / pval)
            track && push!(E.comb, Dict{Int,T}(j => w / pval for (j, w) in comb))
            push!(combsize, size_k / abs(pval))
            push!(E.independent, k)
            row_to_basis[prow] = length(E.pivots)
        end
    end
    return E
end

"""
    PresolveResult{T}

Result of [`SDPAFamily.find_redundant_columns`](@ref).

* `redundant::Vector{Int}`: the sorted indices of the columns that are linear
  combinations of other columns. The corresponding variables are fixed to zero.
* `ray::Union{Nothing,Vector{T}}`: `nothing`, or a vector `d` with `A d ≈ 0` and
  `c'd < 0`, which shows that the dual problem is infeasible.
"""
struct PresolveResult{T}
    redundant::Vector{Int}
    ray::Union{Nothing,Vector{T}}
end

"""
    find_redundant_columns(A::SparseMatrixCSC{T}, c::AbstractVector{T}; tol = default_presolve_tol(T)) -> PresolveResult{T}

Finds the columns of `A` that are linear combinations of earlier columns, using
[`SDPAFamily.eliminate`](@ref). If such a combination is inconsistent with `c`, the
problem `min c'x : A x + b ∈ K` is dual infeasible and the result contains a ray.
"""
function find_redundant_columns(A::SparseMatrixCSC{T}, c::AbstractVector{T};
                                tol = default_presolve_tol(T)) where {T}
    E = eliminate(A, c; tol = tol)
    E.inconsistent == 0 && return PresolveResult{T}(E.redundant, nothing)
    # Repeat the elimination, keeping track of the combinations, to construct the ray.
    E = eliminate(A, c; tol = tol, track = true)
    ray = zeros(T, size(A, 2))
    for (j, w) in E.certificate
        ray[j] = w
    end
    return PresolveResult{T}(Int[], ray)
end

"""
    EqualityElimination{T}

Result of [`SDPAFamily.eliminate_equalities`](@ref). For each `q`, the equality
constraints imply

    x[pivots[q]] = rhs[q] - Σ_j rows[q][j] * x[j],

where all `j` are in `free`.
"""
struct EqualityElimination{T}
    pivots::Vector{Int}
    rows::Vector{Dict{Int,T}}
    rhs::Vector{T}
    # row `q` of the reduced system is `Σ_r comb[q][r] * (row r of A_E)`
    comb::Vector{Dict{Int,T}}
    free::Vector{Int}
    # `nothing`, or `y` with `A_E' y ≈ 0` and `b_E' y < 0`
    infeasible::Union{Nothing,Vector{T}}
end

"""
    eliminate_equalities(AE::SparseMatrixCSC{T}, bE::Vector{T}; tol = default_presolve_tol(T)) -> EqualityElimination{T}

Reduces the system `AE * x + bE = 0` to reduced row echelon form by Gaussian
elimination in the number type `T`.
"""
function eliminate_equalities(AE::SparseMatrixCSC{T}, bE::Vector{T};
                              tol = default_presolve_tol(T)) where {T}
    m, n = size(AE)
    # The columns of the transpose are the rows of `AE`.
    E = eliminate(sparse(transpose(AE)), -bE; tol = tol, track = true)
    if E.inconsistent != 0
        # `d` satisfies `AE' d ≈ 0` and `-bE' d < 0`; a certificate needs `bE' y < 0`.
        y = zeros(T, m)
        for (r, w) in E.certificate
            y[r] = -w
        end
        return EqualityElimination{T}(Int[], Dict{Int,T}[], T[], Dict{Int,T}[],
                                      collect(1:n), y)
    end
    Q = length(E.pivots)
    rows = E.basis
    rhs = E.basis_c
    comb = E.comb
    pivot_index = Dict{Int,Int}(E.pivots[q] => q for q in 1:Q)
    # Back substitution: after the iteration for `s`, row `s` is zero in the pivot
    # columns of all other rows.
    for s in Q-1:-1:1
        row = rows[s]
        for j in collect(keys(row))
            q = get(pivot_index, j, 0)
            q > s || continue
            a = row[j]
            for (l, w) in rows[q]
                row[l] = get(row, l, zero(T)) - a * w
            end
            delete!(row, j)
            rhs[s] -= a * rhs[q]
            for (r, w) in comb[q]
                comb[s][r] = get(comb[s], r, zero(T)) - a * w
            end
        end
        droptol = _drop_tol(T) * maximum(abs, values(row); init = zero(T))
        for (l, w) in collect(row)
            if l != E.pivots[s] && abs(w) <= droptol
                delete!(row, l)
            end
        end
    end
    for q in 1:Q
        delete!(rows[q], E.pivots[q])
    end
    free = sort!(collect(setdiff(1:n, E.pivots)))
    return EqualityElimination{T}(E.pivots, rows, rhs, comb, free, nothing)
end

# Values of all variables, given the values of the free variables. With
# `homogeneous = true` the right-hand side is ignored (for rays).
function expand_solution(elim::EqualityElimination{T}, n::Int, xfree::Vector{T};
                         homogeneous::Bool = false) where {T}
    x = zeros(T, n)
    x[elim.free] = xfree
    for q in eachindex(elim.pivots)
        val = homogeneous ? zero(T) : elim.rhs[q]
        for (j, w) in elim.rows[q]
            val -= w * x[j]
        end
        x[elim.pivots[q]] = val
    end
    return x
end

# Solves `AE' y = g` for `g` in the row space of `AE`.
function equality_duals(elim::EqualityElimination{T}, m::Int, g::Vector{T}) where {T}
    y = zeros(T, m)
    for q in eachindex(elim.pivots)
        β = g[elim.pivots[q]]
        iszero(β) && continue
        for (r, w) in elim.comb[q]
            y[r] += β * w
        end
    end
    return y
end

const SDPAModel{T} = MOI.FileFormats.SDPA.Model{T}
const CONE_TYPES = (MOI.Nonnegatives, MOI.PositiveSemidefiniteConeTriangle)

_cone_constraints(model::SDPAModel{T}, S) where {T} =
    MOI.get(model, MOI.ListOfConstraintIndices{AFF{T},S}())

# The constraint matrix `A` and constant vector `b` of the cone constraints (one row
# per scalar row, in the order in which `MOI.FileFormats.SDPA` writes the blocks, and
# one column per variable), and the objective vector `c` of the minimization problem
# that is written to the SDPA file.
function constraint_matrix(model::SDPAModel{T}) where {T}
    n = MOI.get(model, MOI.NumberOfVariables())
    I, J, V = Int[], Int[], T[]
    b = T[]
    for S in CONE_TYPES
        for ci in _cone_constraints(model, S)
            f = MOI.get(model, MOI.ConstraintFunction(), ci)
            offset = length(b)
            for t in f.terms
                push!(I, offset + t.output_index)
                push!(J, t.scalar_term.variable.value)
                push!(V, t.scalar_term.coefficient)
            end
            append!(b, f.constants)
        end
    end
    A = sparse(I, J, V, length(b), n)
    c = zeros(T, n)
    sense = MOI.get(model, MOI.ObjectiveSense())
    if sense != MOI.FEASIBILITY_SENSE
        obj = MOI.get(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}())
        for t in obj.terms
            c[t.variable.value] += t.coefficient
        end
        if sense == MOI.MAX_SENSE
            c .= .-c
        end
    end
    return A, b, c
end

# A copy of `model` with the variables `keep` (sorted) only. Variable `x[j]` with `j`
# not in `keep` is replaced by `subst[j] = (constant, Dict(i => coefficient))`, an
# affine function of the kept variables, or by zero if `j` is not a key of `subst`.
# The constraints are added in the same order, so the SDPA blocks keep their order.
# Returns the new model and the constant that the substitution adds to the objective.
function substitute_variables(model::SDPAModel{T}, keep::Vector{Int},
                              subst::Dict{Int,Tuple{T,Dict{Int,T}}}) where {T}
    n = MOI.get(model, MOI.NumberOfVariables())
    new_index = zeros(Int, n)
    for (i, j) in enumerate(keep)
        new_index[j] = i
    end
    function substitute(terms)
        constant = zero(T)
        acc = Dict{Int,T}()
        for (j, α) in terms
            if new_index[j] != 0
                acc[new_index[j]] = get(acc, new_index[j], zero(T)) + α
            elseif haskey(subst, j)
                β, expr = subst[j]
                constant += α * β
                for (l, w) in expr
                    acc[new_index[l]] = get(acc, new_index[l], zero(T)) + α * w
                end
            end
        end
        return constant, sort!([(i, v) for (i, v) in acc if !iszero(v)])
    end
    reduced = MOI.FileFormats.SDPA.Model(number_type = T)
    MOI.add_variables(reduced, length(keep))
    for S in CONE_TYPES
        for ci in _cone_constraints(model, S)
            f = MOI.get(model, MOI.ConstraintFunction(), ci)
            constants = copy(f.constants)
            terms = MOI.VectorAffineTerm{T}[]
            rows = [Tuple{Int,T}[] for _ in 1:MOI.output_dimension(f)]
            for t in f.terms
                push!(rows[t.output_index], (t.scalar_term.variable.value, t.scalar_term.coefficient))
            end
            for o in eachindex(rows)
                constant, new_terms = substitute(rows[o])
                constants[o] += constant
                for (i, v) in new_terms
                    push!(terms, MOI.VectorAffineTerm(o,
                        MOI.ScalarAffineTerm(v, MOI.VariableIndex(i))))
                end
            end
            MOI.add_constraint(reduced, MOI.VectorAffineFunction(terms, constants),
                               MOI.get(model, MOI.ConstraintSet(), ci))
        end
    end
    sense = MOI.get(model, MOI.ObjectiveSense())
    MOI.set(reduced, MOI.ObjectiveSense(), sense)
    delta = zero(T)
    if sense != MOI.FEASIBILITY_SENSE
        obj = MOI.get(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}())
        delta, new_terms = substitute([(t.variable.value, t.coefficient) for t in obj.terms])
        MOI.set(reduced, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}(),
                MOI.ScalarAffineFunction(
                    [MOI.ScalarAffineTerm(v, MOI.VariableIndex(i)) for (i, v) in new_terms],
                    zero(T)))
    end
    return reduced, delta
end

"""
    PresolveData{T}

Information that [`SDPAFamily.Optimizer`](@ref) keeps to map the solution of the
problem solved by SDPA back to the model. Variables are numbered as in the
`MOI.FileFormats.SDPA.Model` that is built from the model.
"""
mutable struct PresolveData{T}
    nvars::Int
    # cone constraints and objective of the model, see `constraint_matrix`
    A::SparseMatrixCSC{T,Int}
    b::Vector{T}
    c::Vector{T}
    # equality constraints `AE x + bE = 0`; constraint `k` has rows `eq_rows[k]`
    AE::SparseMatrixCSC{T,Int}
    bE::Vector{T}
    eq_rows::Vector{UnitRange{Int}}
    elim::Union{Nothing,EqualityElimination{T}}
    # variables of the model after the equality elimination that are fixed to zero
    redundant::Vector{Int}
    eq_primal::Vector{T}
    eq_dual::Vector{T}
end

# Block matrices with the given vector (MOI order) in each block.
function blocks_from_vector(blockdims::Vector{Int}, v::Vector{T}) where {T}
    blocks = Matrix{T}[]
    offset = 0
    for d in blockdims
        if d < 0
            B = zeros(T, -d, -d)
            for i in 1:-d
                B[i, i] = v[offset+i]
            end
            offset += -d
        else
            B = zeros(T, d, d)
            k = offset
            for j in 1:d, i in 1:j
                k += 1
                B[i, j] = B[j, i] = v[k]
            end
            offset = k
        end
        push!(blocks, B)
    end
    return blocks
end

# The cone duals as one vector, with the off-diagonal entries of the PSD blocks scaled
# by 2. Then `A' y` is the adjoint of `A` for the inner product `MOI.Utilities.set_dot`
# that MathOptInterface uses for the duals of `PositiveSemidefiniteConeTriangle`.
function weighted_dual_vector(optimizer::Optimizer{T}, M) where {T}
    v = T[]
    for (blk, d) in enumerate(optimizer.blockdims)
        if d < 0
            append!(v, vectorize_block(M, blk, MOI.Nonnegatives))
        else
            B = block(M, blk)
            for j in 1:d, i in 1:j
                push!(v, i == j ? B[i, j] : 2 * B[i, j])
            end
        end
    end
    return v
end

# Whether the constant vector `b` of the cone constraints is in the cone. Returns
# `:feasible` if it is, up to `tol` relative to the largest entry of each block;
# `:nearly` if it is only up to `near_tol`; `:infeasible` otherwise.
function constants_in_cone(blockdims::Vector{Int}, b::Vector{T}, tol, near_tol) where {T}
    result = :feasible
    for B in blocks_from_vector(blockdims, b)
        scale = maximum(abs, B; init = zero(T))
        iszero(scale) && continue
        isposdef(Symmetric(B + tol * scale * I)) && continue
        isposdef(Symmetric(B + near_tol * scale * I)) || return :infeasible
        result = :nearly
    end
    return result
end

# Problem data often comes from `Float64` numbers, so the presolve accepts a point as
# nearly feasible if it is feasible up to the square root of the `Float64` precision.
near_feasibility_tol(::Type{T}) where {T} = max(default_presolve_tol(T), sqrt(eps(Float64)))

"""
    presolve!(optimizer::Optimizer{T}, model::MOI.FileFormats.SDPA.Model{T}, AE, bE, eq_rows) -> Union{Nothing,MOI.FileFormats.SDPA.Model{T}}

Applies the presolve described at the top of `src/presolve.jl`. Returns the reduced
model to pass to SDPA, or `nothing` if the presolve has already determined the result
(an infeasibility certificate, or the solution of a model without remaining
variables); in that case the solution is stored in `optimizer`.
"""
function presolve!(optimizer::Optimizer{T}, model::SDPAModel{T},
                   AE::SparseMatrixCSC{T,Int}, bE::Vector{T},
                   eq_rows::Vector{UnitRange{Int}}) where {T}
    tol = default_presolve_tol(T)
    n = MOI.get(model, MOI.NumberOfVariables())
    A, b, c = constraint_matrix(model)
    data = PresolveData{T}(n, A, b, c, AE, bE, eq_rows, nothing, Int[], T[], T[])
    optimizer.presolve_data = data
    reduced = model
    if !isempty(bE)
        elim = eliminate_equalities(AE, bE; tol = tol)
        data.elim = elim
        if elim.infeasible !== nothing
            set_presolve_result!(optimizer, :presolve_pINF; eq_dual = elim.infeasible)
            return nothing
        end
        subst = Dict{Int,Tuple{T,Dict{Int,T}}}()
        for q in eachindex(elim.pivots)
            subst[elim.pivots[q]] = (elim.rhs[q], Dict{Int,T}(j => -w for (j, w) in elim.rows[q]))
        end
        reduced, delta = substitute_variables(model, elim.free, subst)
        # the objective of the reduced model plus `delta` is the objective of `model`
        optimizer.objconstant += delta
    end
    Ar, br, cr = constraint_matrix(reduced)
    result = find_redundant_columns(Ar, cr; tol = tol)
    if result.ray !== nothing
        x = data.elim === nothing ? result.ray :
            expand_solution(data.elim, n, result.ray; homogeneous = true)
        set_presolve_result!(optimizer, :presolve_dINF; x = x)
        optimizer.primalobj = dot(c, x)
        return nothing
    end
    data.redundant = result.redundant
    nr = size(Ar, 2)
    if length(result.redundant) == nr
        # No variable is left. The only solution is `x = 0` for the reduced model; it is
        # optimal if the constants are in the cone, with the dual solution `y_K = 0`.
        feasible = constants_in_cone(optimizer.blockdims, br, tol, near_feasibility_tol(T))
        phase = feasible == :feasible ? :presolve_OPT :
                feasible == :nearly ? :presolve_NEARLY : :presolve_INF
        set_presolve_result!(optimizer, phase; xr = zeros(T, nr), slack = br)
        return nothing
    end
    if !isempty(result.redundant)
        keep = setdiff(1:nr, result.redundant)
        reduced, _ = substitute_variables(reduced, keep, Dict{Int,Tuple{T,Dict{Int,T}}}())
    end
    return reduced
end

# Stores a result found by the presolve without calling SDPA.
function set_presolve_result!(optimizer::Optimizer{T}, phase::Symbol;
                              x = nothing, xr = nothing, slack = nothing,
                              eq_dual = nothing) where {T}
    data = optimizer.presolve_data
    optimizer.phasevalue = phase
    zero_blocks() = map(n -> zeros(T, abs(n), abs(n)), optimizer.blockdims)
    optimizer.X = PrimalSolution{T}(zero_blocks())
    optimizer.Z = VarDualSolution{T}(slack === nothing ? zero_blocks() :
                                     blocks_from_vector(optimizer.blockdims, slack))
    optimizer.primalobj = zero(T)
    optimizer.dualobj = zero(T)
    if xr !== nothing
        # The objective values of the reduced model are zero; `objconstant` holds the
        # rest of the objective value.
        optimizer.y = xr
        postsolve!(optimizer)
    elseif x !== nothing
        optimizer.y = x
        data.eq_primal = data.AE * x
        data.eq_dual = zeros(T, length(data.bE))
    else
        optimizer.y = zeros(T, data.nvars)
        data.eq_primal = zeros(T, length(data.bE))
        data.eq_dual = eq_dual
        # `MOI.DualObjectiveValue` of the certificate is `-bE' y` for minimization and
        # `bE' y` for maximization, so `objsign * dualobj` with `dualobj = -bE' y`
        optimizer.dualobj = -dot(data.bE, eq_dual)
    end
    return optimizer
end

# Maps the solution of the reduced problem (in `optimizer.y`, with zeros for the
# redundant variables already inserted) to the model.
function postsolve!(optimizer::Optimizer{T}) where {T}
    data = optimizer.presolve_data
    data === nothing && return optimizer
    primal_ray = MOI.get(optimizer, MOI.PrimalStatus()) == MOI.INFEASIBILITY_CERTIFICATE
    dual_ray = MOI.get(optimizer, MOI.DualStatus()) == MOI.INFEASIBILITY_CERTIFICATE
    if data.elim === nothing
        x = optimizer.y
    else
        x = expand_solution(data.elim, data.nvars, optimizer.y; homogeneous = primal_ray)
    end
    optimizer.y = x
    data.eq_primal = primal_ray ? data.AE * x : data.AE * x + data.bE
    if data.elim === nothing
        data.eq_dual = T[]
    else
        yK = weighted_dual_vector(optimizer, optimizer.X)
        g = dual_ray ? -(transpose(data.A) * yK) : data.c - transpose(data.A) * yK
        data.eq_dual = equality_duals(data.elim, length(data.bE), g)
    end
    return optimizer
end
