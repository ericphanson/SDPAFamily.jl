# Part of the presolve (see `src/presolve.jl`): cone membership and the checks of
# certificates, rays and solutions.

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

# Checks a ray `d` of the problem with cone constraints `A x + b ∈ sets`, equality
# constraints `AE x + bE = 0` and objective `c'x`: `AE d = 0`, `A d ∈ K`, `c'd < 0`.
# The checks also reject non-finite values: with `MultiFloat` as the working type, an
# overflow gives `NaN`, for which every comparison is false.
function valid_ray(sets, A, AE, c, d, tol)
    all(isfinite, d) || return false
    rE, mE = _apply(AE, d)
    _all_zero(rE, mE, tol) || return false
    s, m = _apply(A, d)
    in_cone(sets, s, m, tol)[1] || return false
    cd = sum(c .* d; init = zero(eltype(d)))
    return cd < 0 && !_counts_as_zero(cd, sum(abs.(c) .* abs.(d); init = zero(eltype(d))), tol)
end
valid_ray(data::PresolveData, d, tol) = valid_ray(data.sets, data.A, data.AE, data.c, d, tol)

# Checks a certificate of infeasibility (`yK` weighted as in `weighted_dual_vector`,
# `yE`): `A' yK + AE' yE = 0` and `b' yK + bE' yE < 0`. (`yK` is in the dual cone by
# construction for the certificates of the presolve, and for those of SDPA, whose
# dual matrices are positive definite.)
function valid_certificate(A, b, AE, bE, yK, yE, tol)
    all(isfinite, yK) && all(isfinite, yE) || return false
    H = eltype(b)
    r = transpose(A) * yK + transpose(AE) * yE
    m = transpose(abs.(A)) * abs.(yK) + transpose(abs.(AE)) * abs.(yE)
    _all_zero(r, m, tol) || return false
    by = sum(b .* yK; init = zero(H)) + sum(bE .* yE; init = zero(H))
    bm = sum(abs.(b) .* abs.(yK); init = zero(H)) + sum(abs.(bE) .* abs.(yE); init = zero(H))
    return by < 0 && !_counts_as_zero(by, bm, tol)
end
valid_certificate(data::PresolveData, yK, yE, tol) =
    valid_certificate(data.A, data.b, data.AE, data.bE, yK, yE, tol)

# The constraint tolerance `epsilonStar` and the gap tolerance `epsilonDash` of the
# solver parameters, or `nothing` if the parameters are given as a file or as
# `UNSTABLE_BUT_FAST`/`STABLE_BUT_SLOW`.
function solver_tolerances(optimizer::Optimizer{T}) where {T}
    P = optimizer.params
    if P isa ParamsSetting
        P == DEFAULT || return nothing
        P = Params{optimizer.variant,T}()
    end
    return P isa Params ? (Float64(P.epsilonStar), Float64(P.epsilonDash)) : nothing
end
solver_feasibility_tol(optimizer::Optimizer) =
    (t = solver_tolerances(optimizer); t === nothing ? 0.0 : t[1])

# The tolerance of the checks of a solution returned by SDPA against the model: the
# square root of the solver tolerances (or `1e-4` if they are not known). The checks are
# gross: they catch results that are wrong, not results that are inaccurate.
function check_tol(optimizer::Optimizer, data::PresolveData{T,H}) where {T,H}
    t = solver_tolerances(optimizer)
    t === nothing && return max(H(1e-4), sqrt(data.data_tol))
    return sqrt(max(data.data_tol, H(t[1]), H(t[2])))
end

# Gross checks of the result of SDPA, mapped back to the problem with cone constraints
# `A x + b ∈ sets`, equality constraints `AE x + bE = 0` and objective `c'x` (with
# equilibrated rows), at the point `x` with the cone duals `yK` (weighted as in
# `weighted_dual_vector`) and the equality duals `yE`: the primal point (or ray) and
# the duality gap. (Certificates are checked by `check_sdpa_certificate`, and the dual
# equations of the reduced problem in `postsolve!`.) The residuals are compared with
# `tol` times the largest term of the same kind (the rows are equilibrated, so their
# terms are comparable); a residual of an individual row relative to its own terms is
# not meaningful for rows that hold with equality at a point near zero. Returns
# `nothing`, or the failed check.
function check_result(optimizer::Optimizer, sets, A, b, c, AE, bE, x, yK, yE, tol, unit)
    pstat = MOI.get(optimizer, MOI.PrimalStatus())
    dstat = MOI.get(optimizer, MOI.DualStatus())
    H = eltype(b)
    nrm(v) = maximum(abs, v; init = zero(H))
    dot(u, v) = sum(u .* v; init = zero(H))
    all(isfinite, x) && all(isfinite, yK) && all(isfinite, yE) || return "the result is not finite"
    if pstat in (MOI.FEASIBLE_POINT, MOI.INFEASIBILITY_CERTIFICATE)
        # a point, or a ray (`b` and `bE` do not count)
        ray = pstat == MOI.INFEASIBILITY_CERTIFICATE
        s, m = _apply(A, x)
        r, rm = _apply(AE, x)
        if !ray
            s .+= b
            m .+= abs.(b)
            r .+= bE
            rm .+= abs.(bE)
        end
        scale = max(nrm(m), nrm(rm))
        in_cone(sets, s, fill(scale, length(s)), tol)[1] ||
            return "the primal $(ray ? "ray" : "point") violates a cone constraint"
        nrm(r) <= tol * scale || return "the primal $(ray ? "ray" : "point") violates an equality constraint"
        ray && !(dot(c, x) < 0) && return "the primal ray does not decrease the objective"
    end
    if optimizer.phasevalue == :pdOPT
        # The gap relative to its terms, with the floor `unit`: the size in the model of
        # an objective value of 1 in the problem solved by SDPA, for which SDPA's gap
        # criterion is absolute.
        p = dot(c, x)
        d = -(dot(b, yK) + dot(bE, yE))
        mag = dot(abs.(c), abs.(x)) + dot(abs.(b), abs.(yK)) + dot(abs.(bE), abs.(yE))
        abs(p - d) <= tol * (mag + unit) || return "the duality gap is not small"
    end
    return nothing
end

# The tolerance of `check_sdpa_certificate`.
const CERTIFICATE_TOL = 1e-3

# Checks an infeasibility certificate returned by SDPA on the scaled problem written for
# SDPA (`data.sdpa_problem`, before the solution is scaled back), in which solutions of
# size 1 are expected. SDPA returns its last iterate, not a ray, so its certificates
# satisfy `A'y = 0` (or `A d ∈ K`) only approximately; the check asks for the bound
# that the approximate certificate proves:
#
# * dual certificate `y` (`y ∈ K*`, `b'y < 0`, residual `r = A'y`): every feasible `x`
#   has `r'x = y'(A x + b) - b'y ≥ |b'y|`, so `‖x‖∞ ≥ |b'y| / ‖r‖₁`;
# * primal ray `d` (`c'd < 0`, `A d + v ∈ K` with a violation `v`): every dual feasible
#   `y` has `|c'd| = |y'(A d)| ≤ ‖y‖∞ ‖v‖₁` up to the PSD blocks, so `‖y‖∞ ≥ |c'd| / ‖v‖₁`.
#
# The certificate is accepted if this bound is at least `1 / CERTIFICATE_TOL`, that is,
# if it excludes all solutions up to 1000 times the size at which SDPA solved the
# problem. Returns `nothing`, or the failed check.
function check_sdpa_certificate(optimizer::Optimizer{T}) where {T}
    data = optimizer.presolve_data
    pstat = MOI.get(optimizer, MOI.PrimalStatus())
    dstat = MOI.get(optimizer, MOI.DualStatus())
    (pstat == MOI.INFEASIBILITY_CERTIFICATE || dstat == MOI.INFEASIBILITY_CERTIFICATE) || return nothing
    return with_work_precision(T) do
        H = eltype(data.b)
        A, b, c, sets = data.sdpa_problem
        tol = H(CERTIFICATE_TOL)
        dot(u, v) = sum(u .* v; init = zero(H))
        if dstat == MOI.INFEASIBILITY_CERTIFICATE
            y = weighted_dual_vector(optimizer, optimizer.X, H)
            all(isfinite, y) || return "the certificate of infeasibility is not finite"
            by = dot(b, y)
            by < 0 || return "the certificate of infeasibility does not have b'y < 0"
            sum(abs, transpose(A) * y; init = zero(H)) <= tol * abs(by) ||
                return "the certificate of infeasibility does not satisfy A'y = 0 closely enough"
        end
        if pstat == MOI.INFEASIBILITY_CERTIFICATE
            d = H.(optimizer.y[data.keep])
            all(isfinite, d) || return "the primal ray is not finite"
            cd = dot(c, d)
            cd < 0 || return "the primal ray does not decrease the objective"
            # `A d` must be in the cone up to a violation of total size `tol |c'd|`
            s = A * d
            in_cone(sets, s, fill(abs(cd) / length(s), length(s)), tol)[1] ||
                return "the primal ray does not satisfy the cone constraints closely enough"
        end
        return nothing
    end
end
