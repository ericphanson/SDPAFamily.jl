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
# The presolve removes both:
#
# 1. Gaussian elimination on the rows of `A_E` selects one pivot variable per linearly
#    independent row and expresses it as an affine function of the remaining ("free")
#    variables, which is substituted into the cone constraints and the objective.
# 2. Gaussian elimination on the columns of the substituted `A_K` finds a maximal set
#    of linearly independent columns; the other variables are fixed to zero.
#
# ## Arithmetic and decisions
#
# The problem data, given in the number type `T` of the optimizer, are converted exactly
# to `BigFloat` with `2 * precision(T) + 64` bits (the working precision), and the
# presolve computes in that precision. Each decision of the form "this computed number
# is zero" is a one-level check against the data the number was computed from: a
# combination `Σ_j w_j a_j` of data vectors `a_j` with the final coefficients `w_j`
# counts as zero in entry `r` if `|Σ_j w_j a_jr| ≤ τ Σ_j |w_j| |a_jr|`, where
# `τ = sqrt(eps(T) * eps(working precision))`. Rounding errors of the working precision
# are far below `τ`, while numbers that are exact in the data of type `T` (such as an
# exact cancellation `L - L + s` with small `s`) are far above it unless the dynamic
# range within one computed number exceeds about `1/τ`.
#
# * A row of `A_E` (or column of `A_K`) is dependent if all entries of the eliminated
#   vector count as zero; this is confirmed by recomputing the combination directly from
#   the original data. It is consistent if the same combination of the right-hand sides
#   (objective coefficients) counts as zero.
# * After the back substitution, each expression for an eliminated variable is checked
#   against the original rows of `A_E`.
# * A coefficient of the substituted problem is dropped if it counts as zero.
#
# ## Validation and fallback
#
# Every result that the presolve returns without calling SDPA is validated against the
# original data with the same one-level checks: certificates of infeasibility
# (`A'y = 0`, `b'y < 0`, `y` in the dual cone), rays (`A_E d = 0`, `A_K d ∈ K`,
# `c'd < 0`), and the point determined by the equality constraints. After SDPA has
# solved the reduced problem, the solution is mapped back and the original equality
# constraints are checked at the returned point. If any check fails, the presolve falls
# back to passing the problem to SDPA without reduction, with each equality constraint
# written as two inequalities.
#
# ## Mapping the solution back
#
# * The values of the eliminated variables are computed from the free variables.
# * The primal and dual values of the cone constraints are the same in the reduced and
#   the original model.
# * The duals `y_E` of the equality constraints solve `A_E' y_E = c - A_K' y_K` (with
#   the inner product of `MOI.Utilities.set_dot` for the PSD blocks); dropped
#   (dependent) rows get the dual value zero.

const _EXACT = Union{Rational,Integer}

work_type(::Type{T}) where {T} = T <: _EXACT ? T : T === Float64 ? MultiFloats.Float64x3 : BigFloat
work_precision(::Type{T}) where {T} = 2 * precision(T) + 64

# Runs `f()` with the `BigFloat` precision set to the working precision for `T`.
function with_work_precision(f, ::Type{T}) where {T}
    work_type(T) === BigFloat || return f()
    return setprecision(f, BigFloat, work_precision(T))
end

# `eps(T)` in the working precision; `pT` is the precision of `T = BigFloat` of the
# caller (`eps(BigFloat)` depends on the current precision). Call inside
# `with_work_precision`.
eps_of(::Type{T}, pT) where {T} =
    T <: _EXACT ? zero(T) : work_type(T)(2)^(1 - (T === BigFloat ? pT : precision(T)))

# ## Tolerances
#
# The presolve uses two tolerances for its decisions, both relative to the terms that
# are combined (see the description above):
#
# * `rounding_tol` (`τ = sqrt(eps(T) * eps(W))` for the working type `W`): a computed
#   number that is at most `τ` times its terms consists of rounding errors of the
#   working precision. Such numbers are removed during the eliminations and the
#   substitution, and the reduced row echelon form is checked with `τ`.
# * `data_tol` (`64 eps(D)`, where `D` is `Float64` if all data of the model are
#   `Float64` numbers and `T` otherwise): a row or column is dependent, and a dependent
#   row or column is consistent, if the combination is zero up to `data_tol`. Results
#   that rest on these decisions (certificates, rays, the determined point) are
#   validated against the scaled original data with the same tolerance.
#
# The determined point is accepted as feasible with `max(data_tol, epsilonStar)`, the
# constraint tolerance of the solver. After the solve, the residuals of the mapped-back
# solution are compared with the tolerances above and with the residuals of the
# reduced problem (see `postsolve!`).
#
# Call these inside `with_work_precision`.
rounding_tol(::Type{T}, pT) where {T} =
    T <: _EXACT ? zero(T) : sqrt(eps_of(T, pT) * eps(work_type(T)))
data_tol(::Type{T}, pT) where {T} = T <: _EXACT ? zero(T) : 64 * eps_of(T, pT)

_counts_as_zero(value, bound, tol) = abs(value) <= tol * bound

# Converts `x` from the working type to `T`; for `T = BigFloat` with the precision
# `pT` of the caller.
_to(::Type{T}, pT, x) where {T} = T === BigFloat ? BigFloat(x; precision = pT) : T(x)

"""
    Elimination{T}

Result of [`SDPAFamily.eliminate`](@ref).
"""
mutable struct Elimination{T}
    # basis vector `q` has the entry `one(T)` in row `pivots[q]`, and zeros in the pivot
    # rows of the basis vectors before it; it is `Σ_j comb[q][j] * A[:, j]` up to
    # entries that count as zero
    pivots::Vector{Int}
    basis::Vector{Dict{Int,T}}
    basis_c::Vector{T}
    comb::Vector{Dict{Int,T}}
    independent::Vector{Int}
    redundant::Vector{Int}
    # for each redundant column `k`, coefficients `w` (with `w[k] = 1`) such that
    # `A w` is zero up to the dependence tolerance
    dependence::Dict{Int,Dict{Int,T}}
    # the first column that is dependent but inconsistent with `c`, or 0
    inconsistent::Int
    # if `inconsistent > 0`: coefficients `d` with `A d ≈ 0` and `c'd < 0`
    certificate::Dict{Int,T}
    # `true` if a computed result did not pass its check against the data
    failed::Bool
end

# `Σ_j |comb_j| |A[r, j]|` for each row `r` in `rows`
function _bound(A::SparseMatrixCSC{T}, comb::Dict{Int,T}, rows) where {T}
    bound = Dict{Int,T}(r => zero(T) for r in rows)
    rv = rowvals(A)
    nz = nonzeros(A)
    for (j, w) in comb
        aw = abs(w)
        for p in nzrange(A, j)
            r = rv[p]
            haskey(bound, r) && (bound[r] += aw * abs(nz[p]))
        end
    end
    return bound
end

# `Σ_j comb_j A[:, j]` and `Σ_j |comb_j| |A[:, j]|`, computed directly
function _combination(A::SparseMatrixCSC{T}, comb::Dict{Int,T}) where {T}
    val = Dict{Int,T}()
    mag = Dict{Int,T}()
    rv = rowvals(A)
    nz = nonzeros(A)
    for (j, w) in comb
        for p in nzrange(A, j)
            r = rv[p]
            val[r] = get(val, r, zero(T)) + w * nz[p]
            mag[r] = get(mag, r, zero(T)) + abs(w) * abs(nz[p])
        end
    end
    return val, mag
end

"""
    eliminate(A::SparseMatrixCSC{T}, c::AbstractVector{T}; tol = 0, data_tol = tol, pivot_cost = nothing) -> Elimination{T}

Gaussian elimination on the columns of `A`, in the number type `T`. The columns are
processed in order. After the elimination of column `k`, with `w` the coefficients of
the result in terms of the original columns, entry `r` is dropped if its absolute value
is at most `tol * Σ_j |w_j| |A[r, j]| + tol * Σ_j |w_j| max|A[:, j]|` (rounding errors).
Column `k` is linearly dependent on the previous columns if all its entries are at
most this bound with `data_tol` in place of `tol` in the first term, which is confirmed
by recomputing `A w` directly. The same combination of `c` decides whether it is
consistent, with `data_tol`; if not, the elimination stops and records a certificate
`d` with `A d ≈ 0` and `c'd < 0`.

The pivot of a new basis vector is chosen among its entries whose absolute value is at
least `1/10` of the largest one: the entry in the row `r` with the smallest
`pivot_cost[r]` (if given), then the largest entry, then the smallest row.

With `tol = 0` (the default for `Rational` and `Integer` element types) the
computation is exact.
"""
function eliminate(A::SparseMatrixCSC{T}, c::AbstractVector{T};
                   tol = zero(T), data_tol = tol,
                   pivot_cost::Union{Nothing,AbstractVector{Int}} = nothing) where {T}
    n = size(A, 2)
    length(c) == n || throw(DimensionMismatch("`c` must have one entry per column of `A`"))
    colmax = [maximum(abs, view(nonzeros(A), nzrange(A, j)); init = zero(T)) for j in 1:n]
    E = Elimination{T}(Int[], Dict{Int,T}[], T[], Dict{Int,T}[], Int[], Int[],
                       Dict{Int,Dict{Int,T}}(), 0, Dict{Int,T}(), false)
    row_to_basis = Dict{Int,Int}()
    rows_A = rowvals(A)
    vals_A = nonzeros(A)
    for k in 1:n
        v = Dict{Int,T}()
        for p in nzrange(A, k)
            iszero(vals_A[p]) && continue
            r = rows_A[p]
            v[r] = get(v, r, zero(T)) + vals_A[p]
        end
        cv = c[k]
        comb = Dict{Int,T}(k => one(T))
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
                if iszero(newval)
                    delete!(v, r)
                else
                    v[r] = newval
                    q2 = get(row_to_basis, r, 0)
                    q2 > q && push!(todo, q2)
                end
            end
            delete!(v, prow)
            cv -= a * E.basis_c[q]
            for (j, w) in E.comb[q]
                comb[j] = get(comb, j, zero(T)) - a * w
            end
        end
        # entries that count as zero, measured against the original columns
        bound = _bound(A, comb, keys(v))
        # the rounding errors of the working precision in the computation of `v`, which
        # can exceed `tol` times the one-level bound when `comb` has small entries
        spread = tol * sum(abs(w) * colmax[j] for (j, w) in comb)
        for r in collect(keys(v))
            abs(v[r]) <= tol * bound[r] + spread && delete!(v, r)
        end
        if all(abs(v[r]) <= data_tol * bound[r] + spread for r in keys(v))
            # dependent column; confirm with the combination computed directly
            val, mag = _combination(A, comb)
            if !all(abs(val[r]) <= data_tol * mag[r] + spread for r in keys(val))
                E.failed = true
                return E
            end
            cd = zero(T)
            cmag = zero(T)
            for (j, w) in comb
                cd += w * c[j]
                cmag += abs(w) * abs(c[j])
            end
            if abs(cd) > (data_tol + tol) * cmag
                E.inconsistent = k
                s = cd > 0 ? -one(T) : one(T)
                for (j, w) in comb
                    iszero(w) || (E.certificate[j] = s * w)
                end
                return E
            end
            push!(E.redundant, k)
            E.dependence[k] = comb
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
            E.basis[end][prow] = one(T)
            push!(E.basis_c, cv / pval)
            push!(E.comb, Dict{Int,T}(j => w / pval for (j, w) in comb))
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
* `failed::Bool`: `true` if a computed result did not pass its check; then the other
  fields are empty.
"""
struct PresolveResult{T}
    redundant::Vector{Int}
    ray::Union{Nothing,Vector{T}}
    failed::Bool
end

"""
    find_redundant_columns(A::SparseMatrixCSC{T}, c::AbstractVector{T}) -> PresolveResult{T}

Finds the columns of `A` that are linear combinations of earlier columns, using
[`SDPAFamily.eliminate`](@ref). If such a combination is inconsistent with `c`, the
problem `min c'x : A x + b ∈ K` is dual infeasible and the result contains a ray.

The computation is done in the working precision with the tolerances described at the
top of `src/presolve.jl`, and the result is converted to `T`.
"""
function find_redundant_columns(A::SparseMatrixCSC{T}, c::AbstractVector{T}) where {T}
    pT = T === BigFloat ? precision(BigFloat) : 0
    return with_work_precision(T) do
        H = work_type(T)
        r = _find_redundant_columns(H.(A), H.(c), rounding_tol(T, pT), data_tol(T, pT))
        ray = r.ray === nothing ? nothing : T[_to(T, pT, x) for x in r.ray]
        return PresolveResult{T}(r.redundant, ray, r.failed)
    end
end

function _find_redundant_columns(A::SparseMatrixCSC{T}, c::AbstractVector{T}, tol, dtol) where {T}
    E = eliminate(A, c; tol = tol, data_tol = dtol)
    E.failed && return PresolveResult{T}(Int[], nothing, true)
    E.inconsistent == 0 && return PresolveResult{T}(E.redundant, nothing, false)
    ray = zeros(T, size(A, 2))
    for (j, w) in E.certificate
        ray[j] = w
    end
    return PresolveResult{T}(Int[], ray, false)
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
    # the rows of `A_E` that were not found to be dependent, and for each dependent
    # row `r` the coefficients `w` (`w[r] = 1`) of the combination of rows that is zero
    independent::Vector{Int}
    dependence::Dict{Int,Dict{Int,T}}
    # `nothing`, or `y` with `A_E' y ≈ 0` and `b_E' y < 0`
    infeasible::Union{Nothing,Vector{T}}
    # `true` if a computed result did not pass its check against `A_E` and `b_E`
    failed::Bool
end

"""
    eliminate_equalities(AE::SparseMatrixCSC{T}, bE::Vector{T}; pivot_cost = nothing) -> EqualityElimination{T}

Reduces the system `AE * x + bE = 0` to reduced row echelon form by Gaussian
elimination in the number type `T` (see [`SDPAFamily.eliminate`](@ref)), and checks
each resulting row against the original rows. The computation is done in the working
precision and the result is converted to `T`. With `pivot_cost[j]` the number of
nonzeros of variable `j` in the whole model, the pivots are chosen to limit the
fill-in of the substitution.
"""
function eliminate_equalities(AE::SparseMatrixCSC{T}, bE::Vector{T};
                              pivot_cost::Union{Nothing,AbstractVector{Int}} = nothing) where {T}
    pT = T === BigFloat ? precision(BigFloat) : 0
    return with_work_precision(T) do
        H = work_type(T)
        e = _eliminate_equalities(H.(AE), H.(bE), rounding_tol(T, pT), data_tol(T, pT), pivot_cost)
        conv(v) = T[_to(T, pT, x) for x in v]
        convd(d) = Dict{Int,T}(k => _to(T, pT, x) for (k, x) in d)
        return EqualityElimination{T}(e.pivots, convd.(e.rows), conv(e.rhs), convd.(e.comb),
                                      e.free, e.independent,
                                      Dict{Int,Dict{Int,T}}(r => convd(w) for (r, w) in e.dependence),
                                      e.infeasible === nothing ? nothing : conv(e.infeasible),
                                      e.failed)
    end
end

function _eliminate_equalities(AE::SparseMatrixCSC{T}, bE::Vector{T}, tol, dtol,
                               pivot_cost) where {T}
    m, n = size(AE)
    failed() = EqualityElimination{T}(Int[], Dict{Int,T}[], T[], Dict{Int,T}[],
                                      collect(1:n), Int[], Dict{Int,Dict{Int,T}}(), nothing, true)
    # The columns of the transpose are the rows of `AE`.
    AEt = sparse(transpose(AE))
    rowmax = [maximum(abs, view(nonzeros(AEt), nzrange(AEt, r)); init = zero(T)) for r in 1:m]
    E = eliminate(AEt, -bE; tol = tol, data_tol = dtol,
                  pivot_cost = pivot_cost)
    E.failed && return failed()
    if E.inconsistent != 0
        # `d` satisfies `AE' d ≈ 0` and `-bE' d < 0`; a certificate needs `bE' y < 0`.
        y = zeros(T, m)
        for (r, w) in E.certificate
            y[r] = -w
        end
        return EqualityElimination{T}(Int[], Dict{Int,T}[], T[], Dict{Int,T}[],
                                      collect(1:n), Int[], Dict{Int,Dict{Int,T}}(), y, false)
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
                l == j && continue
                newval = get(row, l, zero(T)) - a * w
                iszero(newval) ? delete!(row, l) : (row[l] = newval)
            end
            delete!(row, j)
            rhs[s] -= a * rhs[q]
            for (r, w) in comb[q]
                comb[s][r] = get(comb[s], r, zero(T)) - a * w
            end
        end
    end
    # Drop entries that count as zero, and check each row and right-hand side against
    # the combination of the original rows.
    for q in 1:Q
        val, mag = _combination(AEt, comb[q])
        spread = tol * sum(abs(w) * rowmax[r] for (r, w) in comb[q])
        for l in collect(keys(rows[q]))
            l == E.pivots[q] && continue
            if abs(rows[q][l]) <= tol * get(mag, l, zero(T)) + spread
                delete!(rows[q], l)
            end
        end
        for l in union(keys(val), keys(rows[q]))
            target = l == E.pivots[q] ? one(T) : get(rows[q], l, zero(T))
            if abs(get(val, l, zero(T)) - target) > tol * get(mag, l, zero(T)) + spread
                return failed()
            end
        end
        rd = zero(T)
        rmag = zero(T)
        for (r, w) in comb[q]
            rd -= w * bE[r]
            rmag += abs(w) * abs(bE[r])
        end
        if abs(rd - rhs[q]) > tol * (2 * rmag + sum(abs(w) * rowmax[r] for (r, w) in comb[q]))
            return failed()
        end
        delete!(rows[q], E.pivots[q])
    end
    free = sort!(collect(setdiff(1:n, E.pivots)))
    return EqualityElimination{T}(E.pivots, rows, rhs, comb, free, E.independent,
                                  E.dependence, nothing, false)
end

# Values of all variables, given the values of the free variables. With
# `homogeneous = true` the right-hand side is ignored (for rays).
function expand_solution(elim::EqualityElimination{T}, n::Int, xfree::AbstractVector;
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
function equality_duals(elim::EqualityElimination{T}, m::Int, g::AbstractVector) where {T}
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

# whether all entries are exactly `Float64` numbers
_all_float64(arrays...) = all(all(x -> isfinite(x) && Float64(x) == x, nonzeros_or_all(a)) for a in arrays)
nonzeros_or_all(a::SparseMatrixCSC) = nonzeros(a)
nonzeros_or_all(a) = a

_cone_constraints(model::SDPAModel{T}, S) where {T} =
    MOI.get(model, MOI.ListOfConstraintIndices{AFF{T},S}())

# The constraint matrix `A` and constant vector `b` of the cone constraints (one row
# per scalar row, in the order in which `MOI.FileFormats.SDPA` writes the blocks, and
# one column per variable), the objective vector `c` of the minimization problem that
# is written to the SDPA file, and the list of cone sets in the same order.
function constraint_matrix(model::SDPAModel{T}) where {T}
    n = MOI.get(model, MOI.NumberOfVariables())
    I, J, V = Int[], Int[], T[]
    b = T[]
    sets = Any[]
    for S in CONE_TYPES
        for ci in _cone_constraints(model, S)
            f = MOI.get(model, MOI.ConstraintFunction(), ci)
            push!(sets, MOI.get(model, MOI.ConstraintSet(), ci))
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
    return A, b, c, sets
end

# An SDPA model of type `T` with the cone constraints `A x + b ∈ sets` (rows in the order
# of `constraint_matrix`) and the objective `minimize c'x`. The entries are converted to
# `T` with `_to(T, pT, ⋅)`. Extra `Nonnegatives` constraints `A2 x + b2 ≥ 0` are added
# after the `Nonnegatives` constraints of `sets`.
function build_sdpa_model(::Type{T}, pT, A, b, c, sets; A2 = nothing, b2 = nothing) where {T}
    n = size(A, 2)
    model = MOI.FileFormats.SDPA.Model(number_type = T)
    MOI.add_variables(model, n)
    At = sparse(transpose(A))
    offsets = cumsum(vcat(0, [MOI.dimension(S) for S in sets]))
    function add(Mt, bvec, rows, S)
        terms = MOI.VectorAffineTerm{T}[]
        for (o, r) in enumerate(rows)
            for p in nzrange(Mt, r)
                v = nonzeros(Mt)[p]
                iszero(v) && continue
                push!(terms, MOI.VectorAffineTerm(o, MOI.ScalarAffineTerm(_to(T, pT, v),
                    MOI.VariableIndex(rowvals(Mt)[p]))))
            end
        end
        MOI.add_constraint(model, MOI.VectorAffineFunction(terms, T[_to(T, pT, bvec[r]) for r in rows]), S)
    end
    for (k, S) in enumerate(sets)
        S isa MOI.Nonnegatives && add(At, b, (offsets[k]+1):offsets[k+1], S)
    end
    if A2 !== nothing && !isempty(b2)
        add(sparse(transpose(A2)), b2, 1:length(b2), MOI.Nonnegatives(length(b2)))
    end
    for (k, S) in enumerate(sets)
        S isa MOI.Nonnegatives || add(At, b, (offsets[k]+1):offsets[k+1], S)
    end
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    MOI.set(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}(),
            MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(_to(T, pT, c[j]), MOI.VariableIndex(j))
                                      for j in eachindex(c) if !iszero(c[j])], zero(T)))
    return model
end

# The cone data after substituting the expressions of `elim` for the eliminated
# variables. A coefficient is dropped if it counts as zero against the products it was
# computed from. Returns `A`, `b`, `c` of the reduced problem and the constant `c'x`
# that the substitution adds to the objective.
function substitute(A::SparseMatrixCSC{H}, b::Vector{H}, c::Vector{H},
                    elim::EqualityElimination{H}, tol) where {H}
    n = size(A, 2)
    new_index = zeros(Int, n)
    for (i, j) in enumerate(elim.free)
        new_index[j] = i
    end
    pivot_of = Dict{Int,Int}(p => q for (q, p) in enumerate(elim.pivots))
    val = Dict{Tuple{Int,Int},H}()
    mag = Dict{Tuple{Int,Int},H}()
    bnew = copy(b)
    rv = rowvals(A)
    nz = nonzeros(A)
    for j in 1:n
        for p in nzrange(A, j)
            r, α = rv[p], nz[p]
            if new_index[j] != 0
                key = (r, new_index[j])
                val[key] = get(val, key, zero(H)) + α
                mag[key] = get(mag, key, zero(H)) + abs(α)
            else
                q = pivot_of[j]
                bnew[r] += α * elim.rhs[q]
                for (l, w) in elim.rows[q]
                    key = (r, new_index[l])
                    val[key] = get(val, key, zero(H)) - α * w
                    mag[key] = get(mag, key, zero(H)) + abs(α * w)
                end
            end
        end
    end
    I, J, V = Int[], Int[], H[]
    for (key, v) in val
        (iszero(v) || _counts_as_zero(v, mag[key], tol)) && continue
        push!(I, key[1]); push!(J, key[2]); push!(V, v)
    end
    Anew = sparse(I, J, V, size(A, 1), length(elim.free))
    cnew = c[elim.free]
    cmag = abs.(cnew)
    delta = zero(H)
    for (q, p) in enumerate(elim.pivots)
        iszero(c[p]) && continue
        delta += c[p] * elim.rhs[q]
        for (l, w) in elim.rows[q]
            cnew[new_index[l]] -= c[p] * w
            cmag[new_index[l]] += abs(c[p] * w)
        end
    end
    for i in eachindex(cnew)
        _counts_as_zero(cnew[i], cmag[i], tol) && (cnew[i] = zero(H))
    end
    return Anew, bnew, cnew, delta
end

"""
    Scaling{H}

Power-of-two scaling factors (see [`SDPAFamily.equilibrate`](@ref)). The scaled problem
has the cone constraints `diag(row) * (A * diag(col) * x̂ + rhs * b) ∈ K`, the equality
constraints `diag(eqrow) * (A_E * diag(col) * x̂ + rhs * b_E) = 0` and the objective
`obj * c' * diag(col) * x̂`, with `x = diag(col) * x̂ / rhs`. Since `K` is a cone,
multiplying `b` by `rhs` and the solution by `rhs` gives an equivalent problem.
"""
struct Scaling{H}
    row::Vector{H}
    eqrow::Vector{H}
    col::Vector{H}
    obj::H
    rhs::H
end

# `exponent` of a positive number (for `MultiFloat`, that of its leading limb)
_exponent(x::BigFloat) = exponent(x)
_exponent(x::Rational) = exponent(Float64(x))
_exponent(x) = exponent(Float64(x))
_pow2(::Type{H}, k::Int) where {H} = H(2)^k

# the power of two closest to `1 / sqrt(x)` for `x > 0`
_inv_sqrt_pow2(x::H) where {H} = _pow2(H, -round(Int, _exponent(x) / 2))

# The row group of each cone row: each `Nonnegatives` row is its own group, and the rows
# of a PSD block share one group, so that the scaling of a block is a congruence
# `X ↦ ρ X` and keeps the cone.
function _row_groups(sets, offset::Int)
    groups = Int[]
    g = offset
    for S in sets
        if S isa MOI.Nonnegatives
            for _ in 1:MOI.dimension(S)
                g += 1
                push!(groups, g)
            end
        else
            g += 1
            append!(groups, fill(g, MOI.dimension(S)))
        end
    end
    return groups, g
end

"""
    equilibrate(M::SparseMatrixCSC{H}, groups::Vector{Int}, ngroups::Int; iterations = 20) -> (rowfactors, colfactors)

Ruiz equilibration (D. Ruiz, "A scaling algorithm to equilibrate both rows and columns
norms in matrices", 2001), as used by production LP and conic solvers: the row groups
and the columns of `M` are divided simultaneously by the square roots of their largest
absolute entries, until all of these are within a factor of 2 of 1. The factors are
rounded to powers of two, so scaling by them is exact. Returns the factor of each
group and of each column.
"""
function equilibrate(M::SparseMatrixCSC{H}, groups::Vector{Int}, ngroups::Int;
                     iterations::Int = 20, scale_columns::Bool = true,
                     fixed_column::Union{Nothing,AbstractVector{H}} = nothing) where {H}
    if fixed_column !== nothing
        # a column (the right-hand side) that counts in the row norms but is not scaled
        M = hcat(M, sparse(fixed_column))
    end
    n = size(M, 2)
    g = ones(H, ngroups)
    s = ones(H, n)
    rv = rowvals(M)
    nz = nonzeros(M)
    for _ in 1:iterations
        gnorm = zeros(H, ngroups)
        cnorm = zeros(H, n)
        has_coefficients = falses(ngroups)
        for j in 1:n, p in nzrange(M, j)
            k = groups[rv[p]]
            v = abs(nz[p]) * g[k] * s[j]
            gnorm[k] = max(gnorm[k], v)
            cnorm[j] = max(cnorm[j], v)
            if fixed_column === nothing || j < n
                has_coefficients[k] |= !iszero(nz[p])
            end
        end
        # a row without coefficients is a constant; it is not scaled
        for k in 1:ngroups
            has_coefficients[k] || (gnorm[k] = zero(H))
        end
        balanced(x) = iszero(x) || (x >= one(H) / 2 && x <= 2)
        ncols = fixed_column === nothing ? n : n - 1
        all(balanced, gnorm) && (!scale_columns || all(balanced, view(cnorm, 1:ncols))) && break
        for k in 1:ngroups
            iszero(gnorm[k]) || (g[k] *= _inv_sqrt_pow2(gnorm[k]))
        end
        if scale_columns
            for j in 1:(fixed_column === nothing ? n : n - 1)
                iszero(cnorm[j]) || (s[j] *= _inv_sqrt_pow2(cnorm[j]))
            end
        end
    end
    return g, fixed_column === nothing ? s : s[1:end-1]
end

_scale(M::SparseMatrixCSC, r, c) = Diagonal(r) * M * Diagonal(c)

"""
    PresolveData

Information that `SDPAFamily.Optimizer` keeps to map the solution of the
problem solved by SDPA back to the model. Variables are numbered as in the
`MOI.FileFormats.SDPA.Model` that is built from the model. The arrays of type `H` are in
the working precision.
"""
mutable struct PresolveData{T,H}
    nvars::Int
    pT::Int
    sets::Vector{Any}
    # the scaling before the elimination; the arrays below are scaled by it, and
    # `original` holds the unscaled `(A, b, c, AE, bE)` for the fallback
    pre::Scaling{H}
    original::NTuple{5,Any}
    # the scaling of the reduced problem passed to SDPA, and its columns among the
    # variables of the reduced problem
    post::Union{Nothing,Scaling{H}}
    keep::Vector{Int}
    # the reduced problem (scaled by `pre`, before `post`): cone matrix and objective
    Ared::SparseMatrixCSC{H,Int}
    cred::Vector{H}
    # the dependence and consistency tolerance used (see `data_tol`)
    data_tol::H
    # cone constraints and objective of the scaled model, see `constraint_matrix`
    A::SparseMatrixCSC{H,Int}
    b::Vector{H}
    c::Vector{H}
    # equality constraints `AE x + bE = 0`; constraint `k` has rows `eq_rows[k]`
    AE::SparseMatrixCSC{H,Int}
    bE::Vector{H}
    eq_rows::Vector{UnitRange{Int}}
    elim::Union{Nothing,EqualityElimination{H}}
    # variables of the model after the equality elimination that are fixed to zero
    redundant::Vector{Int}
    # `:reduced` (SDPA solves the reduced problem), `:fallback` (SDPA solves the model
    # with each equality constraint as two inequalities), or `:none` (SDPA not called)
    mode::Symbol
    # the block of the inequalities that replace the equality constraints in `:fallback`
    fallback_block::Int
    fallback_reason::String
    # the constant added to `optimizer.objconstant` by the substitution
    objconstant_delta::T
    eq_primal::Vector{T}
    eq_dual::Vector{T}
end

# Block matrices with the given vector (MOI order) in each block.
function blocks_from_vector(blockdims::Vector{Int}, v::AbstractVector{T}) where {T}
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
function weighted_dual_vector(optimizer::Optimizer, M, ::Type{H}) where {H}
    v = H[]
    for (blk, d) in enumerate(optimizer.blockdims)
        if d < 0
            append!(v, H.(vectorize_block(M, blk, MOI.Nonnegatives)))
        else
            B = block(M, blk)
            for j in 1:d, i in 1:j
                push!(v, i == j ? H(B[i, j]) : 2 * H(B[i, j]))
            end
        end
    end
    return v
end

_block_dims(sets) = [S isa MOI.Nonnegatives ? -MOI.dimension(S) : MOI.side_dimension(S) for S in sets]

"""
    is_psd(B::AbstractMatrix) -> Bool

Whether the symmetric matrix `B` is positive semidefinite, by symmetric Gaussian
elimination with diagonal pivoting. A zero pivot is accepted if the rest of its row is
zero, so rank-deficient matrices (including the zero matrix) are accepted. Exact for
`Rational` element types.
"""
function is_psd(B::AbstractMatrix{T}) where {T}
    M = Matrix{T}(B)
    n = size(M, 1)
    for k in 1:n
        i = k - 1 + argmax([M[l, l] for l in k:n])
        if i != k
            M[[k, i], :] = M[[i, k], :]
            M[:, [k, i]] = M[:, [i, k]]
        end
        d = M[k, k]
        d < 0 && return false
        if iszero(d)
            return all(iszero, @view M[k:n, k:n])
        end
        for r in k+1:n
            l = M[r, k] / d
            iszero(l) && continue
            for s in k+1:n
                M[r, s] -= l * M[k, s]
            end
        end
    end
    return true
end

# Whether `s ∈ K` for the cone blocks `sets`, with magnitudes `m`: a `Nonnegatives`
# entry is accepted if `s[i] ≥ -tol * m[i]`, a PSD block `S` if `S + tol ‖M‖_F I ⪰ 0`.
# Returns `(true, 0)`, `(false, i)` with the most violated `Nonnegatives` entry `i`, or
# `(false, 0)` if only PSD blocks are violated.
function in_cone(sets, s::AbstractVector{H}, m::AbstractVector{H}, tol) where {H}
    worst, worst_ratio = 0, zero(H)
    psd_violated = false
    offset = 0
    for S in sets
        dim = MOI.dimension(S)
        rows = (offset+1):(offset+dim)
        if S isa MOI.Nonnegatives
            for i in rows
                if s[i] < -tol * m[i]
                    ratio = s[i] / m[i]
                    if worst == 0 || ratio < worst_ratio
                        worst, worst_ratio = i, ratio
                    end
                end
            end
        else
            Bs = blocks_from_vector([MOI.side_dimension(S)], s[rows])[1]
            Bm = blocks_from_vector([MOI.side_dimension(S)], m[rows])[1]
            δ = tol * sqrt(sum(abs2, Bm))
            is_psd(Bs + δ * I) || (psd_violated = true)
        end
        offset += dim
    end
    worst > 0 && return false, worst
    return !psd_violated, 0
end

# `M x` and `|M| |x|`
_apply(M, x) = (M * x, abs.(M) * abs.(x))

function _all_zero(val, mag, tol)
    return all(_counts_as_zero(val[i], mag[i], tol) for i in eachindex(val))
end

# Checks a ray `d` against the original data: `AE d = 0`, `A d ∈ K`, `c'd < 0`.
# The checks below also reject non-finite values: with `MultiFloat` as the working type,
# an overflow gives `NaN`, for which every comparison is false.
function valid_ray(data::PresolveData, d, tol)
    all(isfinite, d) || return false
    rE, mE = _apply(data.AE, d)
    _all_zero(rE, mE, tol) || return false
    s, m = _apply(data.A, d)
    in_cone(data.sets, s, m, tol)[1] || return false
    cd = sum(data.c .* d; init = zero(eltype(d)))
    return cd < 0 && !_counts_as_zero(cd, sum(abs.(data.c) .* abs.(d); init = zero(eltype(d))), tol)
end

# Checks a dual certificate (`yK` weighted as in `weighted_dual_vector`, `yE`) against
# the original data: `A' yK + AE' yE = 0`, `yK` in the dual cone, `b' yK + bE' yE < 0`.
function valid_certificate(data::PresolveData, yK, yE, tol)
    all(isfinite, yK) && all(isfinite, yE) || return false
    H = eltype(data.b)
    r = transpose(data.A) * yK + transpose(data.AE) * yE
    m = transpose(abs.(data.A)) * abs.(yK) + transpose(abs.(data.AE)) * abs.(yE)
    _all_zero(r, m, tol) || return false
    by = sum(data.b .* yK; init = zero(H)) + sum(data.bE .* yE; init = zero(H))
    bm = sum(abs.(data.b) .* abs.(yK); init = zero(H)) + sum(abs.(data.bE) .* abs.(yE); init = zero(H))
    return by < 0 && !_counts_as_zero(by, bm, tol)
end

# The constraint tolerance `epsilonStar` of the solver parameters, or zero if the
# parameters are given as a file or as `UNSTABLE_BUT_FAST`/`STABLE_BUT_SLOW`.
function solver_feasibility_tol(optimizer::Optimizer{T}) where {T}
    P = optimizer.params
    if P isa ParamsSetting
        P == DEFAULT || return 0.0
        P = Params{optimizer.variant,T}()
    end
    return P isa Params ? P.epsilonStar : 0.0
end

"""
    presolve!(optimizer::Optimizer{T}, model::MOI.FileFormats.SDPA.Model{T}, AE, bE, eq_rows) -> Union{Nothing,MOI.FileFormats.SDPA.Model{T}}

Applies the presolve described at the top of `src/presolve.jl`. Returns the model to
pass to SDPA (the reduced model, or the fallback model), or `nothing` if the presolve
has already determined a validated result; in that case the result is stored in
`optimizer`.
"""
function presolve!(optimizer::Optimizer{T}, model::SDPAModel{T},
                   AE_T::SparseMatrixCSC{T,Int}, bE_T::Vector{T},
                   eq_rows::Vector{UnitRange{Int}}) where {T}
    pT = T === BigFloat ? precision(BigFloat) : 0
    A_T, b_T, c_T, sets = constraint_matrix(model)
    return with_work_precision(T) do
        H = work_type(T)
        tol = rounding_tol(T, pT)
        ctol = data_tol(T, pT)
        n = size(A_T, 2)
        A0, b0, c0, AE0, bE0 = H.(A_T), H.(b_T), H.(c_T), H.(AE_T), H.(bE_T)
        if !(T <: _EXACT) && T !== Float64 && _all_float64(A_T, b_T, c_T, AE_T, bE_T)
            # data that are exactly `Float64` numbers were most likely computed in
            # `Float64`: decide dependence and consistency at `Float64` precision
            ctol = H(64 * eps(Float64))
        end
        # Scale the rows of the equality and cone constraints (not the variables) before
        # the eliminations, so that the size of a combination of rows is meaningful.
        # Scaling the variables here would let rows that turn out to be redundant
        # distort the scaling of the reduced problem; the variables are scaled after
        # the reduction.
        groups, ng = _row_groups(sets, length(bE0))
        gf, cf = equilibrate(vcat(AE0, A0), vcat(collect(1:length(bE0)), groups), ng;
                             scale_columns = false)
        pre = Scaling{H}(gf[groups], gf[1:length(bE0)], cf, one(H), one(H))
        data = PresolveData{T,H}(n, pT, sets, pre, (A0, b0, c0, AE0, bE0), nothing, Int[],
                                 spzeros(H, 0, 0), H[], ctol,
                                 _scale(A0, pre.row, pre.col), pre.row .* b0, pre.col .* c0,
                                 _scale(AE0, pre.eqrow, pre.col), pre.eqrow .* bE0,
                                 eq_rows, nothing, Int[], :reduced, 0, "", zero(T), T[], T[])
        optimizer.presolve_data = data
        fallback(reason) = fallback_model!(optimizer, reason)
        A, b, c = data.A, data.b, data.c
        delta = zero(H)
        if !isempty(data.bE)
            # the number of nonzeros of each variable, to choose pivots with little fill-in
            pivot_cost = [length(nzrange(A, j)) + length(nzrange(data.AE, j)) for j in 1:n]
            elim = _eliminate_equalities(data.AE, data.bE, tol, ctol, pivot_cost)
            elim.failed && return fallback("the equality elimination did not pass its check")
            if elim.infeasible !== nothing
                y = elim.infeasible
                valid_certificate(data, zeros(H, length(b)), y, ctol) ||
                    return fallback("the certificate for inconsistent equality constraints did not pass its check")
                set_presolve_result!(optimizer, :presolve_pINF; eq_dual = y)
                return nothing
            end
            data.elim = elim
            A, b, c, delta = substitute(A, b, c, elim, tol)
        end
        result = _find_redundant_columns(A, c, tol, ctol)
        redundant = result.failed ? Int[] : result.redundant
        if result.ray !== nothing
            d = data.elim === nothing ? result.ray :
                expand_solution(data.elim, n, result.ray; homogeneous = true)
            valid_ray(data, d, ctol) || return fallback("the ray of the presolve did not pass its check")
            set_presolve_result!(optimizer, :presolve_dINF; x = d)
            optimizer.primalobj = _to(T, pT, sum(data.c .* d; init = zero(H)))
            return nothing
        end
        data.redundant = redundant
        nr = size(A, 2)
        if length(redundant) == nr
            return determined_point!(optimizer, data, delta, tol, ctol)
        end
        keep = setdiff(1:nr, redundant)
        data.objconstant_delta = optimizer.objsign * _to(T, pT, delta)
        optimizer.objconstant += data.objconstant_delta
        # Equilibrate the reduced problem (rows, with the right-hand side counted in
        # the row norms, and columns), then scale the right-hand side and the objective
        # so that their largest entries are in [1, 2), as conic solvers do.
        Ak = A[:, keep]
        groups, ng = _row_groups(sets, 0)
        gf, cf = equilibrate(Ak, groups, ng; fixed_column = b)
        rows = gf[groups]
        bk = rows .* b
        ck = cf .* c[keep]
        to_unit(v) = (m = maximum(abs, v; init = zero(H)); iszero(m) ? one(H) : _pow2(H, -_exponent(m)))
        β = to_unit(bk)
        σ = to_unit(ck)
        data.post = Scaling{H}(rows, H[], cf, σ, β)
        data.keep = keep
        data.Ared = Ak
        data.cred = c[keep]
        return build_sdpa_model(T, pT, _scale(Ak, rows, cf), β .* bk, σ .* ck, sets)
    end
end

# All variables are determined: the equality constraints give the point `x0`, and the
# cone constraints do not depend on the remaining variables. Checks `x0` against the
# original data.
function determined_point!(optimizer::Optimizer{T}, data::PresolveData{T,H}, delta, tol, ctol) where {T,H}
    n = data.nvars
    x0 = data.elim === nothing ? zeros(H, n) :
         expand_solution(data.elim, n, zeros(H, length(data.elim.free)))
    rE, mE = _apply(data.AE, x0)
    _all_zero(rE .+ data.bE, mE .+ abs.(data.bE), ctol) ||
        return fallback_model!(optimizer, "the point determined by the equality constraints did not pass its check")
    s, m = _apply(data.A, x0)
    s .+= data.b
    m .+= abs.(data.b)
    all(isfinite, s) && all(isfinite, m) ||
        return fallback_model!(optimizer, "the determined point is not finite")
    # the objective value at `x0` (the reduced model contributes zero)
    data.objconstant_delta = optimizer.objsign * _to(T, data.pT, delta)
    optimizer.objconstant += data.objconstant_delta
    ftol = max(ctol, H(solver_feasibility_tol(optimizer)))
    feasible, row = in_cone(data.sets, s, m, ftol)
    if feasible
        # optimal with the dual solution `y_K = 0`; the equality duals solve `AE' y_E = c`
        yE = data.elim === nothing ? H[] : equality_duals(data.elim, length(data.bE), data.c)
        r = transpose(data.AE) * yE .- data.c
        rm = transpose(abs.(data.AE)) * abs.(yE) .+ abs.(data.c)
        _all_zero(r, rm, ctol) ||
            return fallback_model!(optimizer, "the dual solution at the determined point did not pass its check")
        set_presolve_result!(optimizer, :presolve_OPT; x = x0, slack = s, eq_dual = yE)
        return nothing
    end
    if row > 0
        # the unit vector of the violated entry, mapped through the equality constraints
        yK = zeros(H, length(s))
        yK[row] = one(H)
        yE = data.elim === nothing ? H[] :
             equality_duals(data.elim, length(data.bE), -(transpose(data.A) * yK))
        if valid_certificate(data, yK, yE, ctol)
            set_presolve_result!(optimizer, :presolve_INF_cert; x = x0, slack = s,
                                 eq_dual = yE, dual_K = yK)
            return nothing
        end
    end
    set_presolve_result!(optimizer, :presolve_INF; x = x0, slack = s)
    return nothing
end

# The cone values of the blocks `M` as one vector in the order of `constraint_matrix`.
function cone_vector(optimizer::Optimizer, M, ::Type{H}) where {H}
    v = H[]
    for (blk, d) in enumerate(optimizer.blockdims)
        S = d < 0 ? MOI.Nonnegatives : MOI.PositiveSemidefiniteConeTriangle
        append!(v, H.(vectorize_block(M, blk, S)))
    end
    return v
end

# Multiplies the cone duals by `dual` and the cone slacks by `slack` (vectors in the
# order of `constraint_matrix`; constant within PSD blocks).
function scale_cone_values!(optimizer::Optimizer{T}, dual, slack, pT) where {T}
    H = eltype(dual)
    conv(v) = T[_to(T, pT, vi) for vi in v]
    optimizer.X = PrimalSolution{T}(blocks_from_vector(optimizer.blockdims,
        conv(cone_vector(optimizer, optimizer.X, H) .* dual)))
    optimizer.Z = VarDualSolution{T}(blocks_from_vector(optimizer.blockdims,
        conv(cone_vector(optimizer, optimizer.Z, H) .* slack)))
    return optimizer
end

# Maps a solution of the scaled model (`data.pre`) to the model.
function unscale_pre!(optimizer::Optimizer{T}) where {T}
    data = optimizer.presolve_data
    with_work_precision(T) do
        H = eltype(data.b)
        conv(v) = T[_to(T, data.pT, vi) for vi in v]
        pre = data.pre
        optimizer.y = conv(pre.col .* H.(optimizer.y))
        scale_cone_values!(optimizer, pre.row, one(H) ./ pre.row, data.pT)
        isempty(data.eq_dual) || (data.eq_dual = conv(pre.eqrow .* H.(data.eq_dual)))
        isempty(data.eq_primal) || (data.eq_primal = conv(H.(data.eq_primal) ./ pre.eqrow))
    end
    return optimizer
end

# Maps the solution of the reduced problem written for SDPA (`data.post`) to the
# variables of the reduced problem.
function unscale_post!(optimizer::Optimizer{T}) where {T}
    data = optimizer.presolve_data
    with_work_precision(T) do
        H = eltype(data.b)
        post = data.post
        y = H.(optimizer.y)
        y[data.keep] .*= post.col ./ post.rhs
        optimizer.y = T[_to(T, data.pT, v) for v in y]
        scale_cone_values!(optimizer, post.row ./ post.obj, one(H) ./ (post.row .* post.rhs), data.pT)
        optimizer.primalobj = _to(T, data.pT, H(optimizer.primalobj) / (post.obj * post.rhs))
        optimizer.dualobj = _to(T, data.pT, H(optimizer.dualobj) / (post.obj * post.rhs))
    end
    return optimizer
end

# Stores a result found by the presolve without calling SDPA. All arguments are in the
# working precision.
function set_presolve_result!(optimizer::Optimizer{T}, phase::Symbol;
                              x = nothing, slack = nothing, eq_dual = nothing,
                              dual_K = nothing) where {T}
    data = optimizer.presolve_data
    H = eltype(data.b)
    conv(v) = T[_to(T, data.pT, vi) for vi in v]
    data.mode = :none
    optimizer.phasevalue = phase
    optimizer.blockdims = _block_dims(data.sets)
    zero_blocks() = map(n -> zeros(T, abs(n), abs(n)), optimizer.blockdims)
    optimizer.X = PrimalSolution{T}(dual_K === nothing ? zero_blocks() :
                                    blocks_from_vector(optimizer.blockdims, conv(dual_K)))
    optimizer.Z = VarDualSolution{T}(slack === nothing ? zero_blocks() :
                                     blocks_from_vector(optimizer.blockdims, conv(slack)))
    optimizer.primalobj = zero(T)
    optimizer.dualobj = zero(T)
    xH = x === nothing ? zeros(H, data.nvars) : x
    optimizer.y = conv(xH)
    if phase == :presolve_dINF
        data.eq_primal = conv(data.AE * xH)
    else
        data.eq_primal = conv(data.AE * xH .+ data.bE)
    end
    yE = eq_dual === nothing ? zeros(H, length(data.bE)) : eq_dual
    data.eq_dual = conv(yE)
    if phase == :presolve_pINF || phase == :presolve_INF_cert
        # `MOI.DualObjectiveValue` of a certificate is `-b'y` for minimization and `b'y`
        # for maximization, so `objsign * dualobj` with `dualobj = -b'y`
        yK = dual_K === nothing ? zeros(H, length(data.b)) : dual_K
        optimizer.dualobj = _to(T, data.pT, -(sum(data.b .* yK; init = zero(H)) +
                                             sum(data.bE .* yE; init = zero(H))))
    end
    # For `:presolve_OPT`, the objective values of the reduced model are zero and
    # `objconstant` holds the objective value.
    return unscale_pre!(optimizer)
end

# Switches to the fallback: the model without reduction, with each equality constraint
# written as two `Nonnegatives` rows, and returns it.
# The number of fallbacks since the package was loaded (the tests check it).
const FALLBACK_COUNT = Ref(0)

function fallback_model!(optimizer::Optimizer{T}, reason::String) where {T}
    FALLBACK_COUNT[] += 1
    data = optimizer.presolve_data
    data.mode = :fallback
    data.fallback_reason = reason
    # the fallback model has the original objective
    optimizer.objconstant -= data.objconstant_delta
    data.objconstant_delta = zero(T)
    data.elim = nothing
    data.redundant = Int[]
    if optimizer.verbosity != SILENT
        @warn "SDPAFamily presolve: $reason. The problem is passed to the solver without the presolve."
    end
    nnn = count(S -> S isa MOI.Nonnegatives, data.sets)
    data.fallback_block = isempty(data.bE) ? 0 : nnn + 1
    A, b, c, AE, bE = data.original
    return with_work_precision(T) do
        build_sdpa_model(T, data.pT, A, b, c, data.sets;
                         A2 = vcat(AE, -AE), b2 = vcat(bE, -bE))
    end
end

# Maps the solution of the problem solved by SDPA (in `optimizer.y`, with zeros for the
# redundant variables already inserted) to the model. In `:reduced` mode, returns
# `false` if the expanded point does not satisfy the original equality constraints.
function postsolve!(optimizer::Optimizer{T}) where {T}
    data = optimizer.presolve_data
    data === nothing && return true
    primal_ray = MOI.get(optimizer, MOI.PrimalStatus()) == MOI.INFEASIBILITY_CERTIFICATE
    dual_ray = MOI.get(optimizer, MOI.DualStatus()) == MOI.INFEASIBILITY_CERTIFICATE
    data.mode == :reduced && unscale_post!(optimizer)
    ok = with_work_precision(T) do
        H = eltype(data.b)
        conv(v) = T[_to(T, data.pT, vi) for vi in v]
        m = length(data.bE)
        if data.mode == :fallback
            _, _, _, AE, bE = data.original
            x = H.(optimizer.y)
            data.eq_primal = conv(primal_ray ? AE * x : AE * x .+ bE)
            if m > 0
                yy = H.(vectorize_block(optimizer.X, data.fallback_block, MOI.Nonnegatives))
                data.eq_dual = conv(yy[1:m] .- yy[m+1:2m])
            else
                data.eq_dual = T[]
            end
            return true
        end
        if data.elim === nothing
            x = H.(optimizer.y)
        else
            x = expand_solution(data.elim, data.nvars, H.(optimizer.y); homogeneous = primal_ray)
        end
        xT = conv(x)
        # Check the original equality constraints at the returned point. The rows that
        # were kept hold up to the rounding of `x` to `T`. A row `i` dropped as
        # dependent holds up to the dependence tolerance times the terms of its
        # combination `w` with other rows, plus the residuals of those rows.
        xr = H.(xT)
        r, rm = _apply(data.AE, xr)
        if !primal_ray
            r .+= data.bE
            rm .+= abs.(data.bE)
        end
        all(isfinite, r) || return false
        ptol = 64 * eps_of(T, data.pT)
        dependence = data.elim === nothing ? Dict{Int,Dict{Int,H}}() : data.elim.dependence
        for i in eachindex(r)
            bound = ptol * rm[i]
            if haskey(dependence, i)
                for (k, w) in dependence[i]
                    k == i && continue
                    bound += abs(w) * (2 * data.data_tol * rm[k] + abs(r[k]))
                end
                bound += 2 * data.data_tol * rm[i]
            end
            abs(r[i]) <= bound || return false
        end
        optimizer.y = xT
        data.eq_primal = conv(r)
        if data.elim === nothing
            data.eq_dual = T[]
        else
            yK = weighted_dual_vector(optimizer, optimizer.X, H)
            g = dual_ray ? -(transpose(data.A) * yK) : data.c - transpose(data.A) * yK
            yE = equality_duals(data.elim, m, g)
            if MOI.get(optimizer, MOI.DualStatus()) == MOI.FEASIBLE_POINT
                # The relative residual of `A' y = c` on the original data must not
                # exceed that of the reduced problem solved by SDPA by much: the mapping
                # back is exact up to rounding, and the dropped columns are consistent
                # up to the dependence tolerance. This catches ill-conditioned
                # eliminations that amplify the solver's error.
                relres(res, scale) = maximum(abs, res; init = zero(H)) /
                                     max(maximum(scale; init = zero(H)), floatmin(Float64))
                res = g .- transpose(data.AE) * yE
                scale = abs.(data.c) .+ transpose(abs.(data.A)) * abs.(yK) .+
                        transpose(abs.(data.AE)) * abs.(yE)
                res_red = data.cred .- transpose(data.Ared) * yK
                scale_red = abs.(data.cred) .+ transpose(abs.(data.Ared)) * abs.(yK)
                relres(res, scale) <= 100 * (relres(res_red, scale_red) + data.data_tol) ||
                    return false
            end
            data.eq_dual = conv(yE)
        end
        return true
    end
    ok && data.mode == :reduced && unscale_pre!(optimizer)
    return ok
end
