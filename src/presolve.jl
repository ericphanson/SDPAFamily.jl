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

# Each number that the presolve computes is a sum of products of the problem data. The
# presolve keeps, for each such number, the sum of the absolute values of these products
# (its "magnitude"). A computed number counts as zero if its absolute value is at most
# `tol` times its magnitude: then it may consist of rounding errors only. Linear
# dependence and consistency are decided with this rule, so a decision does not depend
# on the scale of unrelated entries of the problem.
default_presolve_tol(::Type{T}) where {T} = T <: Union{Rational,Integer} ? zero(T) : 64 * eps(T)

_is_rounding_error(value, magnitude, tol) = abs(value) <= tol * magnitude

"""
    Elimination{T}

Result of [`SDPAFamily.eliminate`](@ref).
"""
mutable struct Elimination{T}
    # basis vector `q` has the entry `one(T)` in row `pivots[q]`, and zeros in the pivot
    # rows of the basis vectors before it; `basis_mag[q]` holds the magnitudes of its
    # entries, and `basis_cmag[q]` the magnitude of `basis_c[q]`
    pivots::Vector{Int}
    basis::Vector{Dict{Int,T}}
    basis_mag::Vector{Dict{Int,T}}
    basis_c::Vector{T}
    basis_cmag::Vector{T}
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
    eliminate(A::SparseMatrixCSC{T}, c::AbstractVector{T}; tol = default_presolve_tol(T), track = false, pivot_cost = nothing) -> Elimination{T}

Gaussian elimination on the columns of `A`, in the number type `T`. The columns are
processed in order. Each entry computed during the elimination counts as zero if its
absolute value is at most `tol` times the sum of the absolute values of the terms
that were added to compute it. Column `k` is linearly dependent on the previous
columns if all its entries become zero. For such a column the same combination is
applied to `c`; if the result does not count as zero by the same rule, column `k` is
recorded as inconsistent and the elimination stops.

The pivot of a new basis vector is chosen among its entries whose absolute value is at
least `1/10` of the largest one: the entry in the row `r` with the smallest
`pivot_cost[r]` (if given), then the largest entry, then the smallest row.

With `tol = 0` (the default for `Rational` and `Integer` element types) the
computation is exact.
"""
function eliminate(A::SparseMatrixCSC{T}, c::AbstractVector{T};
                   tol = default_presolve_tol(T), track::Bool = false,
                   pivot_cost::Union{Nothing,AbstractVector{Int}} = nothing) where {T}
    n = size(A, 2)
    length(c) == n || throw(DimensionMismatch("`c` must have one entry per column of `A`"))
    E = Elimination{T}(Int[], Dict{Int,T}[], Dict{Int,T}[], T[], T[], Dict{Int,T}[],
                       Int[], Int[], 0, Dict{Int,T}())
    row_to_basis = Dict{Int,Int}()
    rows_A = rowvals(A)
    vals_A = nonzeros(A)
    for k in 1:n
        v = Dict{Int,T}()
        vmag = Dict{Int,T}()
        for p in nzrange(A, k)
            iszero(vals_A[p]) && continue
            r = rows_A[p]
            v[r] = get(v, r, zero(T)) + vals_A[p]
            vmag[r] = get(vmag, r, zero(T)) + abs(vals_A[p])
        end
        cv = c[k]
        cmag = abs(c[k])
        comb = track ? Dict{Int,T}(k => one(T)) : Dict{Int,T}()
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
            bmag = E.basis_mag[q]
            for (r, val) in E.basis[q]
                r == prow && continue
                newval = get(v, r, zero(T)) - a * val
                newmag = get(vmag, r, zero(T)) + abs(a) * bmag[r]
                if _is_rounding_error(newval, newmag, tol)
                    delete!(v, r)
                    delete!(vmag, r)
                else
                    v[r] = newval
                    vmag[r] = newmag
                    q2 = get(row_to_basis, r, 0)
                    q2 > q && push!(todo, q2)
                end
            end
            delete!(v, prow)
            delete!(vmag, prow)
            cv -= a * E.basis_c[q]
            cmag += abs(a) * E.basis_cmag[q]
            if track
                for (j, w) in E.comb[q]
                    comb[j] = get(comb, j, zero(T)) - a * w
                end
            end
        end
        if isempty(v)
            # dependent column (this includes a column of zeros)
            if !_is_rounding_error(cv, cmag, tol)
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
            maxabs = maximum(abs, values(v))
            prow, pval = 0, zero(T)
            for (r, val) in v
                abs(val) >= maxabs / 10 || continue
                if prow == 0
                    prow, pval = r, val
                    continue
                end
                if pivot_cost !== nothing && pivot_cost[r] != pivot_cost[prow]
                    pivot_cost[r] < pivot_cost[prow] && ((prow, pval) = (r, val))
                elseif abs(val) > abs(pval) || (abs(val) == abs(pval) && r < prow)
                    prow, pval = r, val
                end
            end
            push!(E.pivots, prow)
            push!(E.basis, Dict{Int,T}(r => val / pval for (r, val) in v))
            push!(E.basis_mag, Dict{Int,T}(r => m / abs(pval) for (r, m) in vmag))
            E.basis[end][prow] = one(T)
            push!(E.basis_c, cv / pval)
            push!(E.basis_cmag, cmag / abs(pval))
            track && push!(E.comb, Dict{Int,T}(j => w / pval for (j, w) in comb))
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

where all `j` are in `free`. `rows_mag` and `rhs_mag` hold the magnitudes of these
numbers (see `default_presolve_tol`).
"""
struct EqualityElimination{T}
    pivots::Vector{Int}
    rows::Vector{Dict{Int,T}}
    rows_mag::Vector{Dict{Int,T}}
    rhs::Vector{T}
    rhs_mag::Vector{T}
    # row `q` of the reduced system is `Σ_r comb[q][r] * (row r of A_E)`
    comb::Vector{Dict{Int,T}}
    free::Vector{Int}
    # `nothing`, or `y` with `A_E' y ≈ 0` and `b_E' y < 0`
    infeasible::Union{Nothing,Vector{T}}
end

"""
    eliminate_equalities(AE::SparseMatrixCSC{T}, bE::Vector{T}; tol = default_presolve_tol(T), pivot_cost = nothing) -> EqualityElimination{T}

Reduces the system `AE * x + bE = 0` to reduced row echelon form by Gaussian
elimination in the number type `T` (see [`SDPAFamily.eliminate`](@ref)). With
`pivot_cost[j]` the number of nonzeros of variable `j` in the whole model, the pivots
are chosen to limit the fill-in of the substitution.
"""
function eliminate_equalities(AE::SparseMatrixCSC{T}, bE::Vector{T};
                              tol = default_presolve_tol(T),
                              pivot_cost::Union{Nothing,AbstractVector{Int}} = nothing) where {T}
    m, n = size(AE)
    # The columns of the transpose are the rows of `AE`.
    E = eliminate(sparse(transpose(AE)), -bE; tol = tol, track = true,
                  pivot_cost = pivot_cost)
    if E.inconsistent != 0
        # `d` satisfies `AE' d ≈ 0` and `-bE' d < 0`; a certificate needs `bE' y < 0`.
        y = zeros(T, m)
        for (r, w) in E.certificate
            y[r] = -w
        end
        return EqualityElimination{T}(Int[], Dict{Int,T}[], Dict{Int,T}[], T[], T[],
                                      Dict{Int,T}[], collect(1:n), y)
    end
    Q = length(E.pivots)
    rows = E.basis
    mags = E.basis_mag
    rhs = E.basis_c
    rhs_mag = E.basis_cmag
    comb = E.comb
    pivot_index = Dict{Int,Int}(E.pivots[q] => q for q in 1:Q)
    # Back substitution: after the iteration for `s`, row `s` is zero in the pivot
    # columns of all other rows.
    for s in Q-1:-1:1
        row = rows[s]
        mag = mags[s]
        for j in collect(keys(row))
            q = get(pivot_index, j, 0)
            q > s || continue
            a = row[j]
            for (l, w) in rows[q]
                l == j && continue
                newval = get(row, l, zero(T)) - a * w
                newmag = get(mag, l, zero(T)) + abs(a) * mags[q][l]
                if _is_rounding_error(newval, newmag, tol)
                    delete!(row, l)
                    delete!(mag, l)
                else
                    row[l] = newval
                    mag[l] = newmag
                end
            end
            delete!(row, j)
            delete!(mag, j)
            rhs[s] -= a * rhs[q]
            rhs_mag[s] += abs(a) * rhs_mag[q]
            for (r, w) in comb[q]
                comb[s][r] = get(comb[s], r, zero(T)) - a * w
            end
        end
    end
    for q in 1:Q
        delete!(rows[q], E.pivots[q])
        delete!(mags[q], E.pivots[q])
    end
    free = sort!(collect(setdiff(1:n, E.pivots)))
    return EqualityElimination{T}(E.pivots, rows, mags, rhs, rhs_mag, comb, free, nothing)
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

# An affine expression `constant + Σ coefficients[i] * x[i]` with the magnitudes of its
# numbers (see `default_presolve_tol`).
struct AffineExpression{T}
    constant::T
    constant_mag::T
    coefficients::Dict{Int,T}
    coefficient_mags::Dict{Int,T}
end

# A copy of `model` with the variables `keep` (sorted) only. Variable `x[j]` with `j`
# not in `keep` is replaced by the affine expression `subst[j]` in the variables of
# `model`, all of which must be in `keep`, or by zero if `j` is not a key of `subst`.
# A coefficient of the result counts as zero if it may consist of rounding errors
# only. The constraints are added in the same order, so the SDPA blocks keep their
# order. Returns the new model, the constant that the substitution adds to the
# objective, and the magnitudes of the constants of the cone constraints (in the order
# of `constraint_matrix`).
function substitute_variables(model::SDPAModel{T}, keep::Vector{Int},
                              subst::Dict{Int,AffineExpression{T}};
                              tol = default_presolve_tol(T)) where {T}
    n = MOI.get(model, MOI.NumberOfVariables())
    new_index = zeros(Int, n)
    for (i, j) in enumerate(keep)
        new_index[j] = i
    end
    function substitute(terms)
        constant = zero(T)
        constant_mag = zero(T)
        acc = Dict{Int,T}()
        accmag = Dict{Int,T}()
        for (j, α) in terms
            if new_index[j] != 0
                i = new_index[j]
                acc[i] = get(acc, i, zero(T)) + α
                accmag[i] = get(accmag, i, zero(T)) + abs(α)
            elseif haskey(subst, j)
                e = subst[j]
                constant += α * e.constant
                constant_mag += abs(α) * e.constant_mag
                for (l, w) in e.coefficients
                    i = new_index[l]
                    acc[i] = get(acc, i, zero(T)) + α * w
                    accmag[i] = get(accmag, i, zero(T)) + abs(α) * e.coefficient_mags[l]
                end
            end
        end
        new_terms = sort!([(i, v) for (i, v) in acc
                           if !iszero(v) && !_is_rounding_error(v, accmag[i], tol)])
        return constant, constant_mag, new_terms
    end
    reduced = MOI.FileFormats.SDPA.Model(number_type = T)
    MOI.add_variables(reduced, length(keep))
    bmag = T[]
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
                constant, constant_mag, new_terms = substitute(rows[o])
                push!(bmag, abs(constants[o]) + constant_mag)
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
        delta, _, new_terms = substitute([(t.variable.value, t.coefficient) for t in obj.terms])
        MOI.set(reduced, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}(),
                MOI.ScalarAffineFunction(
                    [MOI.ScalarAffineTerm(v, MOI.VariableIndex(i)) for (i, v) in new_terms],
                    zero(T)))
    end
    return reduced, delta, bmag
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

# The constraint tolerance `epsilonStar` of the solver parameters, or zero if the
# parameters are given as a file or as `UNSTABLE_BUT_FAST`/`STABLE_BUT_SLOW`.
function solver_feasibility_tol(optimizer::Optimizer{T}) where {T}
    P = optimizer.params
    if P isa ParamsSetting
        P == DEFAULT || return zero(T)
        P = Params{optimizer.variant,T}()
    end
    return P isa Params ? T(P.epsilonStar) : zero(T)
end

# Whether the constant vector `b` of the cone constraints (with magnitudes `bmag`) is
# in the cone. A `Nonnegatives` entry `b[i]` is accepted if `b[i] ≥ -ftol * bmag[i]`.
# A PSD block `B` is accepted if `B + δ I` is positive definite, where `δ` is `ftol`
# times the Frobenius norm of the magnitudes of the block. Returns `(true, 0)`,
# `(false, i)` with the most violated `Nonnegatives` entry `i`, or `(false, 0)` if only
# PSD blocks are violated.
function constants_in_cone(blockdims::Vector{Int}, b::Vector{T}, bmag::Vector{T}, ftol) where {T}
    worst, worst_ratio = 0, zero(T)
    psd_violated = false
    offset = 0
    for d in blockdims
        if d < 0
            for i in offset .+ (1:-d)
                if b[i] < -ftol * bmag[i]
                    ratio = b[i] / bmag[i]
                    if worst == 0 || ratio < worst_ratio
                        worst, worst_ratio = i, ratio
                    end
                end
            end
            offset += -d
        else
            B = zeros(T, d, d)
            M = zeros(T, d, d)
            k = offset
            for j in 1:d, i in 1:j
                k += 1
                B[i, j] = B[j, i] = b[k]
                M[i, j] = M[j, i] = bmag[k]
            end
            offset = k
            δ = ftol * sqrt(sum(abs2, M))
            isposdef(Symmetric(B + δ * I)) || (psd_violated = true)
        end
    end
    worst > 0 && return false, worst
    return !psd_violated, 0
end

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
    bmag = abs.(b)
    if !isempty(bE)
        # the number of nonzeros of each variable, to choose pivots with little fill-in
        pivot_cost = [length(nzrange(A, j)) + length(nzrange(AE, j)) for j in 1:n]
        elim = eliminate_equalities(AE, bE; tol = tol, pivot_cost = pivot_cost)
        data.elim = elim
        if elim.infeasible !== nothing
            set_presolve_result!(optimizer, :presolve_pINF; eq_dual = elim.infeasible)
            return nothing
        end
        subst = Dict{Int,AffineExpression{T}}()
        for q in eachindex(elim.pivots)
            subst[elim.pivots[q]] = AffineExpression{T}(
                elim.rhs[q], elim.rhs_mag[q],
                Dict{Int,T}(j => -w for (j, w) in elim.rows[q]), elim.rows_mag[q])
        end
        reduced, delta, bmag = substitute_variables(model, elim.free, subst; tol = tol)
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
        ftol = max(tol, solver_feasibility_tol(optimizer))
        feasible, row = constants_in_cone(optimizer.blockdims, br, bmag, ftol)
        phase = feasible ? :presolve_OPT : row > 0 ? :presolve_INF_cert : :presolve_INF
        set_presolve_result!(optimizer, phase; xr = zeros(T, nr), slack = br, cert_row = row)
        return nothing
    end
    if !isempty(result.redundant)
        keep = setdiff(1:nr, result.redundant)
        reduced, _, _ = substitute_variables(reduced, keep, Dict{Int,AffineExpression{T}}())
    end
    return reduced
end

# Stores a result found by the presolve without calling SDPA.
function set_presolve_result!(optimizer::Optimizer{T}, phase::Symbol;
                              x = nothing, xr = nothing, slack = nothing,
                              eq_dual = nothing, cert_row::Int = 0) where {T}
    data = optimizer.presolve_data
    optimizer.phasevalue = phase
    zero_blocks() = map(n -> zeros(T, abs(n), abs(n)), optimizer.blockdims)
    optimizer.X = PrimalSolution{T}(zero_blocks())
    optimizer.Z = VarDualSolution{T}(slack === nothing ? zero_blocks() :
                                     blocks_from_vector(optimizer.blockdims, slack))
    optimizer.primalobj = zero(T)
    optimizer.dualobj = zero(T)
    if xr !== nothing
        if cert_row > 0
            # The dual certificate `y_K = e_row` (a violated `Nonnegatives` entry).
            e = zeros(T, length(slack))
            e[cert_row] = one(T)
            optimizer.X = PrimalSolution{T}(blocks_from_vector(optimizer.blockdims, e))
        end
        # The objective values of the reduced model are zero; `objconstant` holds the
        # rest of the objective value.
        optimizer.y = xr
        postsolve!(optimizer)
        if cert_row > 0
            yK = weighted_dual_vector(optimizer, optimizer.X)
            optimizer.dualobj = -(dot(data.b, yK) + dot(data.bE, data.eq_dual))
        end
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
