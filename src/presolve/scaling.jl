# Part of the presolve (see `src/presolve.jl`): the cone data of the SDPA model,
# the substitution, and the scaling.

const SDPAModel{T} = MOI.FileFormats.SDPA.Model{T}
const CONE_TYPES = (MOI.Nonnegatives, MOI.PositiveSemidefiniteConeTriangle)

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
# variables: column `l` of the reduced matrix is the column of the free variable `l`
# minus the columns of the pivot variables times their coefficients of `l`. A
# coefficient that is rounding residue (see `Tolerances`) is dropped. Returns `A`,
# `b`, `c` of the reduced problem and the constant `c'x` that the substitution adds to
# the objective, and the sums `cmag` of the absolute values of the terms of `c`.
function substitute(A::SparseMatrixCSC{H}, b::Vector{H}, c::Vector{H},
                    elim::EqualityElimination{H}, tol::Tolerances{H}) where {H}
    m, n = size(A)
    nfree = length(elim.free)
    new_index = zeros(Int, n)
    for (i, j) in enumerate(elim.free)
        new_index[j] = i
    end
    # for each free variable, the pivot rows `q` in which it occurs, and its coefficients
    occurrences = [Tuple{Int,H}[] for _ in 1:nfree]
    for q in eachindex(elim.pivots), (j, w) in _pairs(elim.rows[q])
        push!(occurrences[new_index[j]], (q, w))
    end
    rv = rowvals(A)
    nz = nonzeros(A)
    acc = Accumulator{H}(m, tol.exact)
    colptr = ones(Int, nfree + 1)
    I, V = Int[], H[]
    cnew = Vector{H}(undef, nfree)
    cmag = Vector{H}(undef, nfree)
    for (i, j) in enumerate(elim.free)
        _reset!(acc)
        cj = c[j]
        cm = abs(cj)
        cemax = iszero(cj) || tol.exact ? 0 : _exponent(cj)
        for p in nzrange(A, j)
            iszero(nz[p]) || _add!(acc, rv[p], nz[p])
        end
        for (q, w) in occurrences[i]
            pv = elim.pivots[q]
            for p in nzrange(A, pv)
                iszero(nz[p]) || _sub!(acc, rv[p], nz[p] * w)
            end
            if !iszero(c[pv])
                t = c[pv] * w
                tol.exact || (cemax = iszero(cj) ? _exponent(t) : max(cemax, _exponent(t)))
                cj -= t
                cm += abs(t)
            end
        end
        for r in sort!(acc.pattern)
            _is_residue(tol, acc.val[r], acc.emax[r]) && continue
            push!(I, r)
            push!(V, acc.val[r])
        end
        colptr[i+1] = length(I) + 1
        cnew[i] = _is_residue(tol, cj, cemax) ? zero(H) : cj
        cmag[i] = cm
    end
    Anew = SparseMatrixCSC(m, nfree, colptr, I, V)
    bnew = copy(b)
    delta = zero(H)
    for (q, pv) in enumerate(elim.pivots)
        iszero(elim.rhs[q]) && continue
        for p in nzrange(A, pv)
            bnew[rv[p]] += nz[p] * elim.rhs[q]
        end
        delta += c[pv] * elim.rhs[q]
    end
    return Anew, bnew, cnew, cmag, delta
end

"""
    Scaling{H}

Power-of-two scaling factors (see [`SDPAFamily.scaled_sdpa_model`](@ref)). The scaled
problem has the cone constraints `diag(row) * (A * diag(col) * x̂ + b) ∈ K` and the
objective `obj * c' * diag(col) * x̂`, with `x = diag(col) * x̂`.
"""
struct Scaling{H}
    row::Vector{H}
    col::Vector{H}
    obj::H
end

_pow2(::Type{H}, k::Int) where {H} = H(2)^k

# `floor(log2(x))` for `x > 0`
_log2(x::BigFloat) = exponent(x)
_log2(x::Rational) = ndigits(numerator(x); base = 2) - ndigits(denominator(x); base = 2)

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
    equilibrate(M::SparseMatrixCSC, groups, ngroups; iterations = 100) -> (groupexponents, colexponents)

Ruiz equilibration (D. Ruiz, "A scaling algorithm to equilibrate both rows and columns
norms in matrices", 2001), the scaling of conic solvers such as SCS, COSMO and
Clarabel: the row groups and the columns of `M` are divided simultaneously by the
square roots of their largest absolute entries, until all of these are in [1/2, 2).
The factors are powers of two, so scaling by them is exact; the function returns their
exponents. Only the binary exponents of the entries are used. Groups and columns
without nonzero entries are not scaled.

The scaling depends only on the largest entries, so entries that are small because
they are rounding errors of the data do not influence it. (Geometric-mean scaling,
which uses the smallest entries too, gave better scalings for sparse matrices whose
variables have very different units, and factors of 2^100 and more for matrices with
such entries.)
"""
function equilibrate(M::SparseMatrixCSC, groups::Vector{Int}, ngroups::Int; iterations::Int = 100)
    n = size(M, 2)
    rv = rowvals(M)
    ex = [iszero(v) ? typemin(Int) : _log2(abs(v)) for v in nonzeros(M)]
    ge = zeros(Int, ngroups)
    se = zeros(Int, n)
    # half of the exponent, rounded away from zero
    half(e) = e > 0 ? cld(e, 2) : -cld(-e, 2)
    for _ in 1:iterations
        gnorm = fill(typemin(Int), ngroups)
        cnorm = fill(typemin(Int), n)
        for j in 1:n, p in nzrange(M, j)
            ex[p] == typemin(Int) && continue
            k = groups[rv[p]]
            v = ex[p] + ge[k] + se[j]
            gnorm[k] = max(gnorm[k], v)
            cnorm[j] = max(cnorm[j], v)
        end
        changed = false
        for (e, norms) in ((ge, gnorm), (se, cnorm)), k in eachindex(e)
            # the largest entry is in [2^x, 2^(x+1))
            x = norms[k]
            (x == typemin(Int) || x == 0 || x == -1) && continue
            e[k] -= half(x)
            changed = true
        end
        changed || break
    end
    return ge, se
end

# The rows of the cone constraints `A x + b ∈ sets` that do not depend on `x`: the rows
# of `Nonnegatives` sets without coefficients, and all rows of PSD blocks without
# coefficients. Returns the rows and, for each, the value of a strictly feasible
# replacement constant (1 for `Nonnegatives` rows, the identity matrix for PSD blocks).
function constant_rows(sets, A::SparseMatrixCSC)
    has_coefficients = falses(size(A, 1))
    rv = rowvals(A)
    nz = nonzeros(A)
    for p in eachindex(rv)
        iszero(nz[p]) || (has_coefficients[rv[p]] = true)
    end
    rows, replacement = Int[], Int[]
    offset = 0
    for S in sets
        dim = MOI.dimension(S)
        if S isa MOI.Nonnegatives
            for i in (offset+1):(offset+dim)
                has_coefficients[i] || (push!(rows, i); push!(replacement, 1))
            end
        elseif !any(view(has_coefficients, (offset+1):(offset+dim)))
            k = offset
            for j in 1:MOI.side_dimension(S), i in 1:j
                k += 1
                push!(rows, k)
                push!(replacement, i == j ? 1 : 0)
            end
        end
        offset += dim
    end
    return rows, replacement
end

# The connected component of each row group and each column in the bipartite graph of
# the nonzero entries of `M` (union-find; the nodes are the groups, then the columns).
function _components(M::SparseMatrixCSC, groups::Vector{Int}, ngroups::Int)
    n = size(M, 2)
    parent = collect(1:(ngroups+n))
    function root(i)
        while parent[i] != i
            parent[i] = parent[parent[i]]
            i = parent[i]
        end
        return i
    end
    rv = rowvals(M)
    nz = nonzeros(M)
    for j in 1:n, p in nzrange(M, j)
        iszero(nz[p]) && continue
        a, b = root(groups[rv[p]]), root(ngroups + j)
        a == b || (parent[a] = b)
    end
    return [root(g) for g in 1:ngroups], [root(ngroups + j) for j in 1:n]
end

"""
    scaled_sdpa_model(T, pT, A, b, c, sets, crows, replacement) -> (model, scaling, (As, bs, cs))

The SDPA model of the problem `minimize c'x` subject to `A x + b ∈ sets`, scaled as
conic solvers such as SCS do:

1. `A` is equilibrated ([`SDPAFamily.equilibrate`](@ref)). The rows of a PSD block
   share one factor, so that the scaling keeps the cone.
2. The right-hand side is normalized. Multiplying the rows of a connected component of
   `A` (of the graph whose edges are its nonzero entries) by `t` and its columns by
   `1/t` leaves `A` unchanged, so the equilibration determines the factors only up to
   one such `t` per component. It is chosen as the power of two that brings the
   largest absolute entry of `b` in the component to [1, 2). For a matrix with one
   component this is one factor for `b`.
3. `c` is multiplied by the power of two that brings its largest absolute entry to
   [1, 2).

The constants of the rows `crows` (rows without coefficients that hold) are replaced
by `replacement`, which keeps the feasible set and makes the rows strictly feasible.
Returns the model, its `Scaling`, and the scaled data `(As, bs, cs)` written to the
model.
"""
function scaled_sdpa_model(::Type{T}, pT, A::SparseMatrixCSC{H}, b::Vector{H}, c::Vector{H},
                           sets, crows, replacement) where {T,H}
    groups, ng = _row_groups(sets, 0)
    ge, se = equilibrate(A, groups, ng)
    # the exponent of the largest entry of the equilibrated `b` in each component
    gcomp, ccomp = _components(A, groups, ng)
    bmax = Dict{Int,Int}()
    for i in eachindex(b)
        iszero(b[i]) && continue
        k = gcomp[groups[i]]
        bmax[k] = max(get(bmax, k, typemin(Int)), _log2(abs(b[i])) + ge[groups[i]])
    end
    for g in 1:ng
        ge[g] -= get(bmax, gcomp[g], 0)
    end
    for j in eachindex(se)
        se[j] += get(bmax, ccomp[j], 0)
    end
    rows = H[_pow2(H, ge[g]) for g in groups]
    cols = H[_pow2(H, e) for e in se]
    cs = cols .* c
    cmax = maximum(abs, cs; init = zero(H))
    σ = iszero(cmax) ? one(H) : _pow2(H, -_log2(cmax))
    bs = rows .* b
    bs[crows] .= H.(replacement)
    As = _scale(A, rows, cols)
    cs .*= σ
    model = build_sdpa_model(T, pT, As, bs, cs, sets)
    return model, Scaling{H}(rows, cols, σ), (As, bs, cs)
end

_scale(M::SparseMatrixCSC, r, c) = Diagonal(r) * M * Diagonal(c)
