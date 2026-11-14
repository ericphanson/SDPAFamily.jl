# Part of the presolve (see `src/presolve.jl`): the working type, the tolerances,
# and the eliminations.

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
# * `data_tol` (`64 eps(T)`, or the `presolve_tol` option of the optimizer if it is
#   larger): a row or column is dependent, and a dependent row or column is consistent,
#   if the combination is zero up to `data_tol`. Results that rest on these decisions
#   (certificates, rays, the determined point) are validated against the scaled
#   original data with the same tolerance. The data are taken as exact numbers of type
#   `T`; for data that were rounded to a coarser type (such as `Float64` numbers in a
#   `BigFloat` model), `presolve_tol` set to `64 eps` of that type lets the presolve
#   treat rows and columns that are dependent up to that rounding as dependent.
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
    # Running error bounds (J. H. Wilkinson; N. J. Higham, "Accuracy and Stability of
    # Numerical Algorithms", 2002, section 3.3): the rounding error of each entry of
    # `basis`, `basis_c` and `comb` is at most a small multiple of `eps(T)` times the
    # corresponding magnitude, which is computed alongside with the absolute values of
    # all operands. A magnitude is at least the absolute value of its entry.
    basis_mag::Vector{Dict{Int,T}}
    basis_c_mag::Vector{T}
    comb_mag::Vector{Dict{Int,T}}
end

# The magnitude of `y - a b`, given magnitudes `ym`, `am`, `bm` of `y`, `a` and `b`
_mag_sub(ym, a, am, b, bm) = ym + abs(a) * bm + am * abs(b)
# The magnitude of `y / p`, given magnitudes `ym` and `pm` of `y` and `p`
_mag_div(y, ym, p, pm) = (ym + abs(y / p) * pm) / abs(p)

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

# `Σ_j comb_j A[:, j]`, `Σ_j |comb_j| |A[:, j]|` and `Σ_j combm_j |A[:, j]|` (a bound
# of the effect of the rounding errors of `comb`), computed directly
function _combination(A::SparseMatrixCSC{T}, comb::Dict{Int,T}, combm::Dict{Int,T}) where {T}
    val = Dict{Int,T}()
    mag = Dict{Int,T}()
    emag = Dict{Int,T}()
    rv = rowvals(A)
    nz = nonzeros(A)
    for (j, w) in comb
        for p in nzrange(A, j)
            r = rv[p]
            val[r] = get(val, r, zero(T)) + w * nz[p]
            mag[r] = get(mag, r, zero(T)) + abs(w) * abs(nz[p])
            emag[r] = get(emag, r, zero(T)) + combm[j] * abs(nz[p])
        end
    end
    return val, mag, emag
end

"""
    eliminate(A::SparseMatrixCSC{T}, c::AbstractVector{T}; tol = 0, data_tol = tol, pivot_cost = nothing) -> Elimination{T}

Gaussian elimination on the columns of `A`, in the number type `T`. The columns are
processed in order. Each computed entry carries a running error bound (a magnitude
`m`, see `Elimination`); an entry `v` with `|v| ≤ tol * m` consists of rounding errors
and is dropped. After the elimination of column `k`, with `w` the coefficients of the
result in terms of the original columns, column `k` is linearly dependent on the
previous columns if each remaining entry satisfies
`|v_r| ≤ data_tol * Σ_j |w_j| |A[r, j]| + tol * m_r`, which is confirmed by recomputing
`A w` directly. The same combination of `c` decides whether it is consistent, with
`data_tol`; if not, the elimination stops and records a certificate `d` with
`A d ≈ 0` and `c'd < 0`.

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
                       Dict{Int,Dict{Int,T}}(), 0, Dict{Int,T}(), false,
                       Dict{Int,T}[], T[], Dict{Int,T}[])
    row_to_basis = Dict{Int,Int}()
    rows_A = rowvals(A)
    vals_A = nonzeros(A)
    for k in 1:n
        # the column, its magnitudes, and the same for `c` and the combination
        v = Dict{Int,T}()
        vm = Dict{Int,T}()
        for p in nzrange(A, k)
            iszero(vals_A[p]) && continue
            r = rows_A[p]
            v[r] = get(v, r, zero(T)) + vals_A[p]
            vm[r] = get(vm, r, zero(T)) + abs(vals_A[p])
        end
        cv, cvm = c[k], abs(c[k])
        comb = Dict{Int,T}(k => one(T))
        combm = Dict{Int,T}(k => one(T))
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
            am = vm[prow]
            bm = E.basis_mag[q]
            for (r, val) in E.basis[q]
                r == prow && continue
                newval = get(v, r, zero(T)) - a * val
                newmag = _mag_sub(get(vm, r, zero(T)), a, am, val, bm[r])
                if iszero(newval)
                    delete!(v, r)
                    delete!(vm, r)
                else
                    v[r] = newval
                    vm[r] = newmag
                    q2 = get(row_to_basis, r, 0)
                    q2 > q && push!(todo, q2)
                end
            end
            delete!(v, prow)
            delete!(vm, prow)
            cvm = _mag_sub(cvm, a, am, E.basis_c[q], E.basis_c_mag[q])
            cv -= a * E.basis_c[q]
            for (j, w) in E.comb[q]
                combm[j] = _mag_sub(get(combm, j, zero(T)), a, am, w, E.comb_mag[q][j])
                comb[j] = get(comb, j, zero(T)) - a * w
            end
        end
        # rounding residue
        for r in collect(keys(v))
            if abs(v[r]) <= tol * vm[r]
                delete!(v, r)
                delete!(vm, r)
            end
        end
        # Dependence. `S = Σ_j |w_j| max|A[:, j]|` bounds the one-level bounds, which
        # are only computed if `v` may be dependent.
        S = sum(abs(w) * colmax[j] for (j, w) in comb)
        dependent = all(abs(v[r]) <= data_tol * S + tol * vm[r] for r in keys(v)) &&
            (bound = _bound(A, comb, keys(v));
             all(abs(v[r]) <= data_tol * bound[r] + tol * vm[r] for r in keys(v)))
        if dependent
            # dependent column; confirm with the combination computed directly
            val, mag, emag = _combination(A, comb, combm)
            if !all(abs(val[r]) <= data_tol * mag[r] + tol * emag[r] for r in keys(val))
                E.failed = true
                return E
            end
            cd = zero(T)
            cmag = zero(T)
            cemag = zero(T)
            for (j, w) in comb
                cd += w * c[j]
                cmag += abs(w) * abs(c[j])
                cemag += combm[j] * abs(c[j])
            end
            if abs(cd) > data_tol * cmag + tol * cemag
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
            pm = vm[prow]
            push!(E.pivots, prow)
            push!(E.basis, Dict{Int,T}(r => val / pval for (r, val) in v))
            push!(E.basis_mag, Dict{Int,T}(r => _mag_div(val, vm[r], pval, pm) for (r, val) in v))
            E.basis[end][prow] = one(T)
            push!(E.basis_c, cv / pval)
            push!(E.basis_c_mag, _mag_div(cv, cvm, pval, pm))
            push!(E.comb, Dict{Int,T}(j => w / pval for (j, w) in comb))
            push!(E.comb_mag, Dict{Int,T}(j => _mag_div(w, combm[j], pval, pm) for (j, w) in comb))
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
    # running error bounds of `rows` and `rhs` (see `Elimination`)
    rows_mag::Vector{Dict{Int,T}}
    rhs_mag::Vector{T}
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
                                      e.failed, convd.(e.rows_mag), conv(e.rhs_mag))
    end
end

function _eliminate_equalities(AE::SparseMatrixCSC{T}, bE::Vector{T}, tol, dtol,
                               pivot_cost) where {T}
    m, n = size(AE)
    failed() = EqualityElimination{T}(Int[], Dict{Int,T}[], T[], Dict{Int,T}[],
                                      collect(1:n), Int[], Dict{Int,Dict{Int,T}}(), nothing, true,
                                      Dict{Int,T}[], T[])
    # The columns of the transpose are the rows of `AE`.
    AEt = sparse(transpose(AE))
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
                                      collect(1:n), Int[], Dict{Int,Dict{Int,T}}(), y, false,
                                      Dict{Int,T}[], T[])
    end
    Q = length(E.pivots)
    rows, rows_mag = E.basis, E.basis_mag
    rhs, rhs_mag = E.basis_c, E.basis_c_mag
    comb, comb_mag = E.comb, E.comb_mag
    pivot_index = Dict{Int,Int}(E.pivots[q] => q for q in 1:Q)
    # Back substitution: after the iteration for `s`, row `s` is zero in the pivot
    # columns of all other rows.
    for s in Q-1:-1:1
        row, rowm = rows[s], rows_mag[s]
        for j in collect(keys(row))
            q = get(pivot_index, j, 0)
            q > s || continue
            a, am = row[j], rowm[j]
            for (l, w) in rows[q]
                l == j && continue
                newval = get(row, l, zero(T)) - a * w
                newmag = _mag_sub(get(rowm, l, zero(T)), a, am, w, rows_mag[q][l])
                if iszero(newval)
                    delete!(row, l)
                    delete!(rowm, l)
                else
                    row[l] = newval
                    rowm[l] = newmag
                end
            end
            delete!(row, j)
            delete!(rowm, j)
            rhs_mag[s] = _mag_sub(rhs_mag[s], a, am, rhs[q], rhs_mag[q])
            rhs[s] -= a * rhs[q]
            for (r, w) in comb[q]
                comb_mag[s][r] = _mag_sub(get(comb_mag[s], r, zero(T)), a, am, w, comb_mag[q][r])
                comb[s][r] = get(comb[s], r, zero(T)) - a * w
            end
        end
    end
    # Drop the rounding residue of the back substitution.
    for q in 1:Q
        for l in collect(keys(rows[q]))
            if l == E.pivots[q] || abs(rows[q][l]) <= tol * rows_mag[q][l]
                delete!(rows[q], l)
                delete!(rows_mag[q], l)
            end
        end
    end
    free = sort!(collect(setdiff(1:n, E.pivots)))
    # A posteriori check of the reduction (a Freivalds-type test, O(nnz) instead of
    # comparing every reduced row with its combination of original rows): at a point
    # with free values spread over [1, 2), the rows that were not found dependent must
    # hold, up to `sqrt(tol * dtol)` relative to their terms. The values are the
    # fractional parts of multiples of the golden ratio, which are reproducible and do
    # not repeat. (The dependent rows were checked when they were found.)
    result = EqualityElimination{T}(E.pivots, rows, rhs, comb, free, E.independent,
                                    E.dependence, nothing, false, rows_mag, rhs_mag)
    phi = (sqrt(5.0) - 1) / 2
    x = expand_solution(result, n, [one(T) + T(mod(k * phi, 1.0)) for k in 1:length(free)])
    r, rm = AE * x .+ bE, abs.(AE) * abs.(x) .+ abs.(bE)
    check_tol = T <: _EXACT ? zero(T) : sqrt(tol * dtol)
    all(abs(r[i]) <= check_tol * rm[i] for i in E.independent) || return failed()
    return result
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
