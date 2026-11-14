# Part of the presolve (see `src/presolve.jl`): the working type, the tolerances,
# and the eliminations.

const _EXACT = Union{Rational,Integer}

# The working type: `BigFloat` with `work_precision` bits, or `T` itself for exact types.
work_type(::Type{T}) where {T} = T <: _EXACT ? T : BigFloat
# Twice the precision of the data, so that sums of products of two data numbers are
# exact, plus guard bits.
work_precision(p::Integer) = 2 * p + 64

# The precision of `T`; `pT` is the precision of `BigFloat` of the caller.
_precision(::Type{T}, pT) where {T} = T === BigFloat ? pT : precision(T)

# Runs `f()` with the `BigFloat` precision set to the working precision for `T`.
function with_work_precision(f, ::Type{T}) where {T}
    T <: _EXACT && return f()
    return setprecision(f, BigFloat, work_precision(precision(T)))
end

# `eps(T)` in the working type. Call inside `with_work_precision`.
eps_of(::Type{T}, pT) where {T} =
    T <: _EXACT ? zero(T) : BigFloat(2)^(1 - _precision(T, pT))

# Converts `x` from the working type to `T`; for `T = BigFloat` with the precision
# `pT` of the caller.
_to(::Type{T}, pT, x) where {T} = T === BigFloat ? BigFloat(x; precision = pT) : T(x)

"""
    Tolerances{H}

The tolerances of the presolve for the working type `H` (see the description at the top
of `src/presolve.jl`).

* `drop_bits`: a computed number whose leading bit is more than `drop_bits` bits below
  the leading bit of the largest term it was computed from (catastrophic cancellation
  of the terms) is rounding residue of the working precision and is set to zero.
* `data`: the tolerance `δ` for linear dependence and consistency, relative to the
  absolute values of the terms of a combination of data vectors.
* `suspect_bits`: a reduced vector all of whose entries lost at least this many bits
  is tested for dependence at `δ`.

For exact types (`exact = true`) only numbers that are exactly zero count as zero.
"""
struct Tolerances{H}
    exact::Bool
    drop_bits::Int
    suspect_bits::Int
    data::H
end

# Call inside `with_work_precision`. `data_tol` is the `presolve_tol` option.
function Tolerances(::Type{T}, pT; data_tol = nothing) where {T}
    T <: _EXACT && return Tolerances{T}(true, 0, 0, zero(T))
    p = _precision(T, pT)
    δ = 64 * eps_of(T, pT)
    data_tol === nothing || (δ = max(δ, BigFloat(data_tol)))
    return Tolerances{BigFloat}(false, (p + work_precision(p)) ÷ 2, max(1, -exponent(δ) ÷ 2), δ)
end

# The binary exponent of a nonzero number (only called for inexact types).
_exponent(x::BigFloat) = exponent(x)
_exponent(x) = exponent(Float64(x))

# Whether the computed number `v`, whose largest term has the exponent `emax`, is
# rounding residue.
_is_residue(tol::Tolerances, v, emax::Int) =
    iszero(v) || (!tol.exact && _exponent(v) < emax - tol.drop_bits)

_counts_as_zero(value, bound, tol) = abs(value) <= tol * bound

# A sparse vector: indices (not sorted) and values.
struct Entries{H}
    idx::Vector{Int}
    val::Vector{H}
end
Entries{H}() where {H} = Entries{H}(Int[], H[])
Base.length(e::Entries) = length(e.idx)
_pairs(e::Entries) = zip(e.idx, e.val)

# A sparse accumulator (J. R. Gilbert, C. Moler, R. Schreiber, "Sparse matrices in
# MATLAB: design and implementation", 1992): a dense vector of values with the list
# of the positions in use. It also records the exponent of the largest term added to
# each position, for `_is_residue`.
struct Accumulator{H}
    val::Vector{H}
    emax::Vector{Int}
    marked::Vector{Bool}
    pattern::Vector{Int}
    exact::Bool
end
Accumulator{H}(n::Int, exact::Bool) where {H} =
    Accumulator{H}(Vector{H}(undef, n), zeros(Int, n), fill(false, n), Int[], exact)

# Adds the nonzero term `t` at position `r`; returns `true` if the position is new.
@inline function _add!(acc::Accumulator, r::Int, t)
    e = acc.exact ? 0 : _exponent(t)
    if acc.marked[r]
        acc.val[r] += t
        e > acc.emax[r] && (acc.emax[r] = e)
        return false
    end
    acc.marked[r] = true
    push!(acc.pattern, r)
    acc.val[r] = t
    acc.emax[r] = e
    return true
end

function _reset!(acc::Accumulator)
    for r in acc.pattern
        acc.marked[r] = false
    end
    empty!(acc.pattern)
    return acc
end

"""
    Elimination{H}

Result of [`SDPAFamily.eliminate`](@ref): an LU factorization of the independent
columns ("basis vectors") of a matrix `A`, in the order in which they were found.

Basis vector `q` is column `vec[q]` of `A`. Its reduced form
`U[q] = A[:, vec[q]] - Σ_{q' < q} L[q][q'] * U[q']` has the entry `pivot[q]` in row
`coord[q]` and is zero in the rows `coord[q']` for `q' < q`; `Uc[q]` is the same
combination of the entries of `c`.
"""
mutable struct Elimination{H}
    vec::Vector{Int}
    coord::Vector{Int}
    pivot::Vector{H}
    U::Vector{Entries{H}}
    Uc::Vector{H}
    L::Vector{Entries{H}}
    # the columns found dependent on earlier columns and consistent with `c`, and for
    # each such column `k` the coefficients `w` (with `w[k] = 1`) such that `A w` is
    # zero up to the dependence tolerance
    redundant::Vector{Int}
    dependence::Dict{Int,Dict{Int,H}}
    # the first column that is dependent but inconsistent with `c`, or 0
    inconsistent::Int
    # if `inconsistent > 0`: coefficients `d` with `A d ≈ 0` and `c'd < 0`
    certificate::Dict{Int,H}
    # `true` if the working precision was not sufficient (two ways of computing the
    # same vector disagree)
    failed::Bool
end

Elimination{H}() where {H} =
    Elimination{H}(Int[], Int[], H[], Entries{H}[], H[], Entries{H}[], Int[],
                   Dict{Int,Dict{Int,H}}(), 0, Dict{Int,H}(), false)

# The threshold of the pivot choice: a pivot is at least this fraction of the largest
# entry of its vector. 0.1 is the usual value of sparse LU codes (MA28, UMFPACK, LUSOL).
const PIVOT_THRESHOLD = 10

"""
    eliminate(A::SparseMatrixCSC{H}, c::AbstractVector{H}, tol::Tolerances{H}; pivot_cost = nothing) -> Elimination{H}

Gaussian elimination on the columns of `A` in the number type `H`, column by column
(a left-looking sparse LU factorization with threshold pivoting, as in J. R. Gilbert
and T. Peierls, "Sparse partial pivoting in time proportional to arithmetic
operations", 1988, with the roles of rows and columns exchanged).

Each column is reduced by the basis vectors found before it. Entries of the result
that are rounding residue (see [`SDPAFamily.Tolerances`](@ref)) are removed. If nothing
remains, or if every remaining entry is small compared to its terms, the column is
tested for dependence on the data: with `w` the coefficients of the reduced column in
terms of the columns of `A` (from a triangular solve), the column is dependent if
`|A w| ≤ δ |A| |w|` holds for every entry, with `δ = tol.data`. This is the
componentwise relative backward error of W. Oettli and W. Prager (1964): some matrix
whose entries differ from those of `A` by relative perturbations of at most `δ` has
`w` in its null space. The test does not depend on the scaling of the rows or columns
of `A`. A dependent column is consistent with `c` if `|c'w| ≤ δ |c|'|w|`; if not, the
elimination stops and records the certificate `d = ±w` with `c'd < 0`.

The pivot of a new basis vector is chosen among its entries whose absolute value is at
least `1/10` of the largest one: the entry in the row `r` with the smallest
`pivot_cost[r]` (if given), then the largest entry, then the smallest row.
"""
function eliminate(A::SparseMatrixCSC{H}, c::AbstractVector{H}, tol::Tolerances{H};
                   pivot_cost::Union{Nothing,AbstractVector{Int}} = nothing) where {H}
    m, n = size(A)
    length(c) == n || throw(DimensionMismatch("`c` must have one entry per column of `A`"))
    E = Elimination{H}()
    coord_to_basis = zeros(Int, m)
    acc = Accumulator{H}(m, tol.exact)
    rows_A = rowvals(A)
    vals_A = nonzeros(A)
    todo = BitSet()
    for k in 1:n
        _reset!(acc)
        empty!(todo)
        for p in nzrange(A, k)
            iszero(vals_A[p]) && continue
            r = rows_A[p]
            _add!(acc, r, vals_A[p])
            q = coord_to_basis[r]
            q > 0 && push!(todo, q)
        end
        cv = c[k]
        cemax = iszero(cv) || tol.exact ? 0 : _exponent(cv)
        l = Entries{H}()
        # Reduce with the basis vectors in the order in which they were found: basis
        # vector `q` is zero in the pivot rows of the basis vectors before it.
        while !isempty(todo)
            q = first(todo)
            delete!(todo, q)
            prow = E.coord[q]
            a = acc.val[prow]
            acc.val[prow] = zero(H)
            _is_residue(tol, a, acc.emax[prow]) && continue
            λ = a / E.pivot[q]
            push!(l.idx, q)
            push!(l.val, λ)
            u = E.U[q]
            for i in eachindex(u.idx)
                r = u.idx[i]
                r == prow && continue
                if _add!(acc, r, -(λ * u.val[i]))
                    q2 = coord_to_basis[r]
                    q2 > 0 && push!(todo, q2)
                end
            end
            if !iszero(E.Uc[q])
                t = λ * E.Uc[q]
                tol.exact || (cemax = iszero(cv) ? _exponent(t) : max(cemax, _exponent(t)))
                cv -= t
            end
        end
        _is_residue(tol, cv, cemax) && (cv = zero(H))
        # the reduced column without rounding residue
        v = Entries{H}()
        suspect = true
        for r in acc.pattern
            x = acc.val[r]
            (coord_to_basis[r] > 0 || _is_residue(tol, x, acc.emax[r])) && continue
            push!(v.idx, r)
            push!(v.val, x)
            suspect = suspect && !tol.exact && _exponent(x) < acc.emax[r] - tol.suspect_bits
        end
        if suspect
            w = _combination(E, l, k, tol)
            if _is_dependent(A, w, tol)
                cd = zero(H)
                cmag = zero(H)
                for (j, wj) in w
                    cd += wj * c[j]
                    cmag += abs(wj) * abs(c[j])
                end
                if !_counts_as_zero(cd, cmag, tol.data)
                    E.inconsistent = k
                    s = cd > 0 ? -one(H) : one(H)
                    for (j, wj) in w
                        E.certificate[j] = s * wj
                    end
                    return E
                end
                push!(E.redundant, k)
                E.dependence[k] = w
                continue
            elseif isempty(v.idx)
                # the reduced column is rounding residue, but the same combination
                # computed from the data is not zero
                E.failed = true
                return E
            end
        end
        prow, pval = _choose_pivot(v, pivot_cost)
        push!(E.vec, k)
        push!(E.coord, prow)
        push!(E.pivot, pval)
        push!(E.U, v)
        push!(E.Uc, cv)
        push!(E.L, l)
        coord_to_basis[prow] = length(E.vec)
    end
    return E
end

function _choose_pivot(v::Entries{H}, pivot_cost) where {H}
    maxabs = maximum(abs, v.val)
    prow, pval = 0, zero(H)
    for (r, val) in _pairs(v)
        PIVOT_THRESHOLD * abs(val) >= maxabs || continue
        if prow == 0
            prow, pval = r, val
        elseif pivot_cost !== nothing && pivot_cost[r] != pivot_cost[prow]
            pivot_cost[r] < pivot_cost[prow] && ((prow, pval) = (r, val))
        elseif abs(val) > abs(pval) || (abs(val) == abs(pval) && r < prow)
            prow, pval = r, val
        end
    end
    return prow, pval
end

# Solves `L' z = l` in place (`z` holds `l` on entry), for the unit lower triangular
# matrix `L` of the elimination. With `tol`, entries of `z` that are rounding residue
# are set to zero.
function _solve_Lt!(z::Vector{H}, E::Elimination{H}, tol::Union{Nothing,Tolerances{H}} = nothing) where {H}
    track = tol !== nothing && !tol.exact
    emax = track ? [iszero(x) ? 0 : _exponent(x) for x in z] : Int[]
    for q in length(z):-1:1
        zq = z[q]
        iszero(zq) && continue
        if track && _is_residue(tol, zq, emax[q])
            z[q] = zero(H)
            continue
        end
        for (q2, λ) in _pairs(E.L[q])
            t = λ * zq
            if track
                e = _exponent(t)
                emax[q2] = iszero(z[q2]) && emax[q2] == 0 ? e : max(emax[q2], e)
            end
            z[q2] -= t
        end
    end
    return z
end

# The coefficients `w` (with `w[k] = 1`) of the column `k` reduced with the multipliers
# `l`, in terms of the columns of `A`: `A[:, k] - Σ_q l[q] U[q] = A w`.
function _combination(E::Elimination{H}, l::Entries{H}, k::Int, tol::Tolerances{H}) where {H}
    z = zeros(H, length(E.vec))
    for (q, λ) in _pairs(l)
        z[q] = λ
    end
    _solve_Lt!(z, E, tol)
    w = Dict{Int,H}(k => one(H))
    for q in eachindex(z)
        iszero(z[q]) || (w[E.vec[q]] = -z[q])
    end
    return w
end

# Whether `|A w| ≤ δ |A| |w|` holds for every entry, computed from the data.
function _is_dependent(A::SparseMatrixCSC{H}, w::Dict{Int,H}, tol::Tolerances{H}) where {H}
    val = Dict{Int,H}()
    mag = Dict{Int,H}()
    rv = rowvals(A)
    nz = nonzeros(A)
    for (j, wj) in w
        for p in nzrange(A, j)
            r = rv[p]
            t = wj * nz[p]
            val[r] = get(val, r, zero(H)) + t
            mag[r] = get(mag, r, zero(H)) + abs(t)
        end
    end
    return all(_counts_as_zero(val[r], mag[r], tol.data) for r in keys(val))
end

"""
    PresolveResult{T}

Result of [`SDPAFamily.find_redundant_columns`](@ref).

* `redundant::Vector{Int}`: the sorted indices of the columns that are linear
  combinations of earlier columns. The corresponding variables are fixed to zero.
* `ray::Union{Nothing,Vector{T}}`: `nothing`, or a vector `d` with `A d ≈ 0` and
  `c'd < 0`, which shows that the dual problem is infeasible.
* `failed::Bool`: `true` if the working precision was not sufficient; then the other
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
        r = _find_redundant_columns(H.(A), H.(c), Tolerances(T, pT))
        ray = r.ray === nothing ? nothing : T[_to(T, pT, x) for x in r.ray]
        return PresolveResult{T}(r.redundant, ray, r.failed)
    end
end

function _find_redundant_columns(A::SparseMatrixCSC{H}, c::AbstractVector{H}, tol::Tolerances{H}) where {H}
    E = eliminate(A, c, tol)
    E.failed && return PresolveResult{H}(Int[], nothing, true)
    E.inconsistent == 0 && return PresolveResult{H}(E.redundant, nothing, false)
    ray = zeros(H, size(A, 2))
    for (j, w) in E.certificate
        ray[j] = w
    end
    return PresolveResult{H}(Int[], ray, false)
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
    rows::Vector{Entries{T}}
    rhs::Vector{T}
    free::Vector{Int}
    # the rows of `A_E` that were not found to be dependent, and for each dependent
    # row `r` the coefficients `w` (`w[r] = 1`) of the combination of rows that is zero
    independent::Vector{Int}
    dependence::Dict{Int,Dict{Int,T}}
    # `nothing`, or `y` with `A_E' y ≈ 0` and `b_E' y < 0`
    infeasible::Union{Nothing,Vector{T}}
    # `true` if the working precision was not sufficient
    failed::Bool
    # the factorization of the independent rows (the columns of `A_E'`), for the duals
    factors::Elimination{T}
end

_failed_elimination(::Type{T}, n; infeasible = nothing, failed = true) where {T} =
    EqualityElimination{T}(Int[], Entries{T}[], T[], collect(1:n), Int[],
                           Dict{Int,Dict{Int,T}}(), infeasible, failed, Elimination{T}())

"""
    eliminate_equalities(AE::SparseMatrixCSC{T}, bE::Vector{T}; pivot_cost = nothing) -> EqualityElimination{T}

Solves the system `AE * x + bE = 0` for one variable per linearly independent row, by
Gaussian elimination on the rows of `AE` (see [`SDPAFamily.eliminate`](@ref)) followed
by back substitution, and checks the result against `AE` and `bE`. The computation is
done in the working precision and the result is converted to `T`. With `pivot_cost[j]`
the number of nonzeros of variable `j` in the whole model, the pivots are chosen to
limit the fill-in of the substitution.
"""
function eliminate_equalities(AE::SparseMatrixCSC{T}, bE::Vector{T};
                              pivot_cost::Union{Nothing,AbstractVector{Int}} = nothing) where {T}
    pT = T === BigFloat ? precision(BigFloat) : 0
    return with_work_precision(T) do
        H = work_type(T)
        e = _eliminate_equalities(H.(AE), H.(bE), Tolerances(T, pT), pivot_cost)
        conv(v) = T[_to(T, pT, x) for x in v]
        conve(v) = Entries{T}(v.idx, conv(v.val))
        convd(d) = Dict{Int,T}(k => _to(T, pT, x) for (k, x) in d)
        F = e.factors
        factors = Elimination{T}(F.vec, F.coord, conv(F.pivot), conve.(F.U), conv(F.Uc), conve.(F.L),
                                 F.redundant, Dict{Int,Dict{Int,T}}(), 0, Dict{Int,T}(), F.failed)
        return EqualityElimination{T}(e.pivots, conve.(e.rows), conv(e.rhs), e.free, e.independent,
                                      Dict{Int,Dict{Int,T}}(r => convd(w) for (r, w) in e.dependence),
                                      e.infeasible === nothing ? nothing : conv(e.infeasible),
                                      e.failed, factors)
    end
end

function _eliminate_equalities(AE::SparseMatrixCSC{H}, bE::Vector{H}, tol::Tolerances{H},
                               pivot_cost) where {H}
    m, n = size(AE)
    # The columns of the transpose are the rows of `AE`; row `k` reads `AE[k, :] x = -bE[k]`.
    AEt = sparse(transpose(AE))
    E = eliminate(AEt, -bE, tol; pivot_cost = pivot_cost)
    E.failed && return _failed_elimination(H, n)
    if E.inconsistent != 0
        # `d` satisfies `AE' d ≈ 0` and `-bE' d < 0`; a certificate needs `bE' y < 0`.
        y = zeros(H, m)
        for (r, w) in E.certificate
            y[r] = -w
        end
        return _failed_elimination(H, n; infeasible = y, failed = false)
    end
    # Back substitution. Reduced row `q` reads `Σ_j U[q][j] x[j] = Uc[q]`, and is zero
    # in the pivot variables of the rows before it. Going backwards, the pivot
    # variables of the later rows are replaced by their expressions.
    Q = length(E.vec)
    is_pivot = zeros(Int, n)
    for q in 1:Q
        is_pivot[E.coord[q]] = q
    end
    rows = Vector{Entries{H}}(undef, Q)
    rhs = Vector{H}(undef, Q)
    acc = Accumulator{H}(n, tol.exact)
    for q in Q:-1:1
        _reset!(acc)
        β = E.Uc[q]
        βemax = iszero(β) || tol.exact ? 0 : _exponent(β)
        for (j, u) in _pairs(E.U[q])
            s = is_pivot[j]
            if s == 0
                _add!(acc, j, u)
            elseif s != q
                for (l, w) in _pairs(rows[s])
                    _add!(acc, l, -(u * w))
                end
                if !iszero(rhs[s])
                    t = u * rhs[s]
                    tol.exact || (βemax = iszero(β) ? _exponent(t) : max(βemax, _exponent(t)))
                    β -= t
                end
            end
        end
        row = Entries{H}()
        for l in acc.pattern
            _is_residue(tol, acc.val[l], acc.emax[l]) && continue
            push!(row.idx, l)
            push!(row.val, acc.val[l] / E.pivot[q])
        end
        rows[q] = row
        rhs[q] = _is_residue(tol, β, βemax) ? zero(H) : β / E.pivot[q]
    end
    free = [j for j in 1:n if is_pivot[j] == 0]
    result = EqualityElimination{H}(E.coord, rows, rhs, free, E.vec, E.dependence, nothing, false, E)
    # A posteriori check of the reduction (a test of R. Freivalds' kind, O(nnz)): at a
    # point with free values spread over [1, 2), the rows that were not found
    # dependent must hold, up to the geometric mean of the rounding and the data
    # tolerance, relative to their terms. The values are the fractional parts of
    # multiples of the golden ratio, which are reproducible and do not repeat. (The
    # dependent rows were checked when they were found.)
    phi = (sqrt(5.0) - 1) / 2
    x = expand_solution(result, n, [one(H) + H(mod(k * phi, 1.0)) for k in 1:length(free)])
    r, rm = AE * x .+ bE, abs.(AE) * abs.(x) .+ abs.(bE)
    check_tol = tol.exact ? zero(H) : sqrt(H(2)^(-tol.drop_bits) * tol.data)
    all(abs(r[i]) <= check_tol * rm[i] for i in E.vec) || return _failed_elimination(H, n)
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
        for (j, w) in _pairs(elim.rows[q])
            val -= w * x[j]
        end
        x[elim.pivots[q]] = val
    end
    return x
end

# Solves `AE' y = g` for `g` in the row space of `AE`, with `y` zero in the dependent
# rows: the independent rows are `L U` (see `Elimination`), so `U' t = g` (forward
# substitution on the pivot variables) and `L' y = t`.
function equality_duals(elim::EqualityElimination{T}, m::Int, g::AbstractVector) where {T}
    E = elim.factors
    Q = length(E.vec)
    res = Vector{T}(g)
    t = zeros(T, Q)
    for q in 1:Q
        tq = res[E.coord[q]] / E.pivot[q]
        iszero(tq) && continue
        t[q] = tq
        for (j, u) in _pairs(E.U[q])
            res[j] -= u * tq
        end
    end
    _solve_Lt!(t, E)
    y = zeros(T, m)
    for q in 1:Q
        y[E.vec[q]] = t[q]
    end
    return y
end
