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

# The tolerances `epsilonStar` and `epsilonDash` of the solver parameters. The preset
# parameter modes of SDPA use `1e-7` for both. If the parameters are given as a file
# that cannot be read as a parameter file of SDPA, returns `nothing`.
function solver_tolerances(optimizer::Optimizer{T}) where {T}
    P = optimizer.params
    if P isa ParamsSetting
        P == DEFAULT || return (1e-7, 1e-7)
        P = Params{optimizer.variant,T}()
    end
    P isa Params && return (Float64(P.epsilonStar), Float64(P.epsilonDash))
    try
        # lines 2 and 10 of a parameter file of SDPA
        lines = readlines(P)
        value(i) = parse(Float64, first(split(lines[i])))
        return (value(2), value(10))
    catch
        return nothing
    end
end
solver_feasibility_tol(optimizer::Optimizer) =
    (t = solver_tolerances(optimizer); t === nothing ? 0.0 : t[1])

# The tolerance of `check_sdpa_result`: the square root of the solver tolerances, that
# is, half of the digits that SDPA was asked for (the convention of `isapprox`); `1e-4`
# if they are not known. The check catches results that are wrong, not results that
# are inaccurate.
function check_tol(optimizer::Optimizer, data::PresolveData{T,H}) where {T,H}
    t = solver_tolerances(optimizer)
    t === nothing && return max(H(1e-4), sqrt(data.data_tol))
    return sqrt(max(data.data_tol, H(t[1]), H(t[2])))
end

# The tolerance of the check of an infeasibility certificate returned by SDPA, see
# `check_sdpa_result`.
const CERTIFICATE_TOL = 1e-3

"""
    check_sdpa_result(optimizer) -> Union{Nothing,String}

Checks the result of SDPA on the scaled problem `min c'x : A x + b ∈ K` that was
written for SDPA (`data.sdpa_problem`), before the solution is scaled back. Returns
`nothing`, or the failed check.

A primal or dual solution that SDPA reports as feasible, and a result that is reported
as optimal, must satisfy the stopping criteria of SDPA with the tolerance `γ` of
`check_tol` in place of `epsilonDash` and `epsilonStar`:

* `A x + b + γ s I ∈ K` with `s = max(1, ‖|A| |x| + |b|‖∞)`,
* `‖c - A'y‖∞ ≤ γ max(1, ‖|c| + |A|'|y|‖∞)`,
* `|c'x + b'y| ≤ γ max(1, (|c'x| + |b'y|) / 2)` (for a result reported as optimal,
  which includes SDPA's status `pdFEAS`).

The gap criterion must also hold in the units of the model, for the gap that SDPA's
tolerance admits: with `p` and `d` the objective values of the scaled problem, `σ` the
factor of the objective, and `p_m`, `d_m` the objective values of the model (`p / σ`,
`d / σ` plus the constant of the substitution),

    max(|p - d|, epsilonStar max(1, |p|, |d|)) / σ ≤ γ max(1, (|p_m| + |d_m|) / 2).

SDPA applies its gap tolerance to the scaled objective, with the absolute floor 1. If
the scaled optimal value is far below 1 and `1/σ` is large (the scaling could not
tell that the objective is small), a gap within the tolerance is large compared to
the objective value of the model; without the presolve SDPA would not have accepted
it.

SDPA does not return rays. It reports unboundedness when the objective of a feasible
iterate passes its parameter `lowerBound` or `upperBound` (`±1e5` by default), and
returns that iterate, which satisfies `A'y = c` (or `A d + b ∈ K`) instead of
`A'y = 0` (`A d ∈ K`). The check is the criterion of the solvers that detect
infeasibility from approximate certificates (SCS, OSQP, Clarabel), with the
tolerance `ε = CERTIFICATE_TOL`:

* dual certificate `y` (`y ∈ K*`): `b'y < 0` and `‖A'y‖∞ ≤ ε |b'y|`. Every feasible
  `x` has `(A'y)'x = y'(A x + b) - b'y ≥ |b'y|`, so `‖x‖₁ ≥ 1/ε`;
* primal ray `d`: `c'd < 0` and `A d + ε |c'd| I ∈ K` (`k ε` for a PSD block of side
  dimension `k`). Every dual feasible `y` has
  `|c'd| = |y'(A d)| ≤ ε |c'd| ‖y‖₁` (with the trace norm for PSD blocks), so
  `‖y‖₁ ≥ 1/ε`.

In the scaled problem `‖c‖∞ < 2` and `‖b‖∞ < 2`, so an iterate that SDPA reports with
its default bounds has `‖A'y‖∞ / |b'y| < 2e-5`. The tolerance `1e-3` accepts
certificates that are 50 times weaker, and rejects those that only exclude solutions
of norm less than 1000 in the scaled problem, in which solutions of norm 1 are
expected.
"""
function check_sdpa_result(optimizer::Optimizer{T}) where {T}
    data = optimizer.presolve_data
    pstat = MOI.get(optimizer, MOI.PrimalStatus())
    dstat = MOI.get(optimizer, MOI.DualStatus())
    return with_work_precision(T) do
        H = eltype(data.b)
        A, b, c, sets = data.sdpa_problem
        γ = check_tol(optimizer, data)
        ε = H(CERTIFICATE_TOL)
        nrm(v) = maximum(abs, v; init = zero(H))
        dot(u, v) = sum(u .* v; init = zero(H))
        x = H.(optimizer.y[data.keep])
        y = weighted_dual_vector(optimizer, optimizer.X, H)
        all(isfinite, x) && all(isfinite, y) || return "the result of SDPA is not finite"
        if pstat == MOI.FEASIBLE_POINT
            s, m = _apply(A, x)
            s .+= b
            m .+= abs.(b)
            in_cone(sets, s, fill(max(one(H), nrm(m)), length(s)), γ)[1] ||
                return "the primal point violates a cone constraint"
        elseif pstat == MOI.INFEASIBILITY_CERTIFICATE
            cd = dot(c, x)
            cd < 0 || return "the primal ray does not decrease the objective"
            # `A d` must be in the cone up to a violation of size `ε |c'd|`
            s = A * x
            in_cone(sets, s, fill(abs(cd), length(s)), ε)[1] ||
                return "the primal ray does not satisfy the cone constraints closely enough"
        end
        if dstat == MOI.FEASIBLE_POINT
            r = c .- transpose(A) * y
            m = abs.(c) .+ transpose(abs.(A)) * abs.(y)
            nrm(r) <= γ * max(one(H), nrm(m)) || return "the dual solution violates the dual equations"
        elseif dstat == MOI.INFEASIBILITY_CERTIFICATE
            by = dot(b, y)
            by < 0 || return "the certificate of infeasibility does not have b'y < 0"
            nrm(transpose(A) * y) <= ε * abs(by) ||
                return "the certificate of infeasibility does not satisfy A'y = 0 closely enough"
        end
        if MOI.get(optimizer, MOI.TerminationStatus()) == MOI.OPTIMAL
            p, d = dot(c, x), -dot(b, y)
            abs(p - d) <= γ * max(one(H), (abs(p) + abs(d)) / 2) || return "the duality gap is not small"
            # The gap in the units of the model. SDPA accepts a gap of `εgap` times
            # `max(1, objective)` of the scaled problem, so the gap can be that large.
            t = solver_tolerances(optimizer)
            εgap = t === nothing ? γ^2 : H(t[1])
            gap = max(abs(p - d), εgap * max(one(H), abs(p), abs(d)))
            σ = data.post.obj
            shift = optimizer.objsign * H(data.objconstant_delta)
            pm, dm = p / σ + shift, d / σ + shift
            gap / σ <= γ * max(one(H), (abs(pm) + abs(dm)) / 2) ||
                return "the gap tolerance of SDPA does not resolve the objective value of the model"
        end
        return nothing
    end
end
