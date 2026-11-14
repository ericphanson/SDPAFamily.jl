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
# to a working type with more precision (`Float64x3` for `Float64`, otherwise `BigFloat`
# with `2 * precision(T) + 64` bits), and the presolve computes in that type. Each
# computed number carries a running error bound (a magnitude `m` computed alongside with
# the absolute values of the operands). A number `v` with `|v| ≤ τ m`, where
# `τ = sqrt(eps(T) * eps(working type))`, consists of rounding errors and is dropped.
# Rounding errors of the working type are far below `τ m`, while numbers that are exact
# in the data of type `T` (such as an exact cancellation `L - L + s` with small `s`, or a
# small coefficient next to a large one) are far above it unless the dynamic range
# within one computation exceeds about `1/τ`.
#
# * A row of `A_E` (or column of `A_K`) is dependent if every entry `v_r` of the
#   eliminated vector satisfies `|v_r| ≤ δ Σ_j |w_j| |a_jr| + τ m_r` for the
#   combination `w` of the data vectors `a_j`, where `δ = data_tol`; this is confirmed
#   by recomputing the combination directly from the data. It is consistent if the
#   same combination of the right-hand sides (objective coefficients) is zero in the
#   same sense.
# * After the back substitution, the reduced system is checked against the original
#   rows of `A_E` at one point.
#
# ## Validation, scaling and fallback
#
# Every result that the presolve returns without calling SDPA is validated against the
# original data at `δ`: certificates of infeasibility (`A'y = 0`, `b'y < 0`), rays
# (`A_E d = 0`, `A_K d ∈ K`, `c'd < 0`), and the point determined by the equality
# constraints. Rows without coefficients are checked as constants and replaced by
# strictly feasible ones. The reduced problem is scaled (Ruiz equilibration of the
# matrix, then one factor each for the right-hand side and the objective). After SDPA
# has solved it, the solution is mapped back, and the original equality constraints,
# the dual residual, and gross conditions on the result (cone constraints, dual
# equations, certificates, duality gap) are checked. If any check fails, the presolve
# falls back to passing the problem to SDPA without reduction, with each equality
# constraint written as two inequalities and the same scaling; a result of the
# fallback that fails the gross checks gets the status `NUMERICAL_ERROR`.
#
# ## Mapping the solution back
#
# * The values of the eliminated variables are computed from the free variables.
# * The primal and dual values of the cone constraints are the same in the reduced and
#   the original model.
# * The duals `y_E` of the equality constraints solve `A_E' y_E = c - A_K' y_K` (with
#   the inner product of `MOI.Utilities.set_dot` for the PSD blocks); dropped
#   (dependent) rows get the dual value zero.

include("presolve/elimination.jl")
include("presolve/scaling.jl")

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
    # the scaling of the problem passed to SDPA, and its columns among the variables of
    # the reduced problem (all variables in `:fallback` mode)
    post::Union{Nothing,Scaling{H}}
    keep::Vector{Int}
    # the rows without coefficients whose constants were replaced (see
    # `scaled_sdpa_model`), and their values before the scaling `post`
    const_rows::Vector{Int}
    const_vals::Vector{H}
    # the scaled data `(A, b, c)` and the cone sets of the problem written for SDPA
    sdpa_problem::Any
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
    # the status reported by SDPA and the failed check, if the result did not pass the
    # checks after the fallback (status `:unverified`)
    sdpa_phase::Symbol
    check_failure::String
end

include("presolve/checks.jl")

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
        if optimizer.presolve_tol !== nothing && !(T <: _EXACT)
            ctol = max(ctol, H(optimizer.presolve_tol))
        end
        n = size(A_T, 2)
        A0, b0, c0, AE0, bE0 = H.(A_T), H.(b_T), H.(c_T), H.(AE_T), H.(bE_T)
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
                                 Int[], H[], nothing, spzeros(H, 0, 0), H[], ctol,
                                 _scale(A0, pre.row, pre.col), pre.row .* b0, pre.col .* c0,
                                 _scale(AE0, pre.eqrow, pre.col), pre.eqrow .* bE0,
                                 eq_rows, nothing, Int[], :reduced, 0, "", zero(T), T[], T[],
                                 :not_called, "")
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
        keep = setdiff(1:size(A, 2), redundant)
        data.objconstant_delta = optimizer.objsign * _to(T, pT, delta)
        optimizer.objconstant += data.objconstant_delta
        if isempty(keep)
            return determined_point!(optimizer, data, tol, ctol)
        end
        Ak = A[:, keep]
        # Rows without coefficients are constants: check them at the point given by the
        # equality constraints (any point, since the rows do not depend on the free
        # variables), and replace them by strictly feasible constants. Otherwise SDPA
        # sees a constraint without interior (for a zero constant) or a constant whose
        # size distorts the scaling.
        crows, replacement = constant_rows(sets, Ak)
        if !isempty(crows)
            x0, s, m = equality_point(data)
            feasible, row = in_cone(sets, _restrict(s, crows), _restrict(m, crows), feasibility_tol(optimizer, data))
            feasible || return infeasible_result!(optimizer, data, x0, s, row)
            data.const_rows = crows
            data.const_vals = s[crows]
        end
        model, post, scaled = scaled_sdpa_model(T, pT, Ak, b, c[keep], sets, crows, replacement)
        data.post = post
        data.sdpa_problem = (scaled..., sets)
        data.keep = keep
        data.Ared = Ak
        data.cred = c[keep]
        return model
    end
end

# `v` with the entries outside `rows` set to zero (a zero block is in every cone).
function _restrict(v::Vector{H}, rows) where {H}
    w = zeros(H, length(v))
    w[rows] .= v[rows]
    return w
end

# The constraint tolerance for decisions on the constants: `max(data_tol, epsilonStar)`.
feasibility_tol(optimizer::Optimizer, data::PresolveData{T,H}) where {T,H} =
    max(data.data_tol, H(solver_feasibility_tol(optimizer)))

# The point `x0` given by the equality constraints with the free variables at zero,
# and the values `s = A x0 + b` and magnitudes `m` of the cone rows.
function equality_point(data::PresolveData{T,H}) where {T,H}
    n = data.nvars
    x0 = data.elim === nothing ? zeros(H, n) :
         expand_solution(data.elim, n, zeros(H, length(data.elim.free)))
    s, m = _apply(data.A, x0)
    s .+= data.b
    m .+= abs.(data.b)
    return x0, s, m
end

# All variables are determined: the equality constraints give the point `x0`, and the
# cone constraints do not depend on the remaining variables. Checks `x0` against the
# original data.
function determined_point!(optimizer::Optimizer{T}, data::PresolveData{T,H}, tol, ctol) where {T,H}
    x0, s, m = equality_point(data)
    rE, mE = _apply(data.AE, x0)
    _all_zero(rE .+ data.bE, mE .+ abs.(data.bE), ctol) ||
        return fallback_model!(optimizer, "the point determined by the equality constraints did not pass its check")
    all(isfinite, s) && all(isfinite, m) ||
        return fallback_model!(optimizer, "the determined point is not finite")
    feasible, row = in_cone(data.sets, s, m, feasibility_tol(optimizer, data))
    feasible || return infeasible_result!(optimizer, data, x0, s, row)
    # optimal with the dual solution `y_K = 0`; the equality duals solve `AE' y_E = c`
    yE = data.elim === nothing ? H[] : equality_duals(data.elim, length(data.bE), data.c)
    r = transpose(data.AE) * yE .- data.c
    rm = transpose(abs.(data.AE)) * abs.(yE) .+ abs.(data.c)
    _all_zero(r, rm, data.data_tol) ||
        return fallback_model!(optimizer, "the dual solution at the determined point did not pass its check")
    set_presolve_result!(optimizer, :presolve_OPT; x = x0, slack = s, eq_dual = yE)
    return nothing
end

# A row of the cone constraints that does not depend on the free variables is violated
# at `x0` (`row > 0` for a `Nonnegatives` row, `row == 0` for a PSD block). For a
# `Nonnegatives` row, the unit vector of the row, mapped through the equality
# constraints, is a certificate of infeasibility.
function infeasible_result!(optimizer::Optimizer{T}, data::PresolveData{T,H}, x0, s, row) where {T,H}
    if row > 0
        yK = zeros(H, length(s))
        yK[row] = one(H)
        yE = data.elim === nothing ? H[] :
             equality_duals(data.elim, length(data.bE), -(transpose(data.A) * yK))
        if valid_certificate(data, yK, yE, data.data_tol)
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
# order of `constraint_matrix`; constant within PSD blocks). With `rows`, sets the
# duals of these rows to zero and their slacks to `values` afterwards.
function scale_cone_values!(optimizer::Optimizer{T}, dual, slack, pT;
                            rows = Int[], values = eltype(dual)[]) where {T}
    H = eltype(dual)
    conv(v) = T[_to(T, pT, vi) for vi in v]
    X = cone_vector(optimizer, optimizer.X, H) .* dual
    Z = cone_vector(optimizer, optimizer.Z, H) .* slack
    X[rows] .= zero(H)
    Z[rows] .= values
    optimizer.X = PrimalSolution{T}(blocks_from_vector(optimizer.blockdims, conv(X)))
    optimizer.Z = VarDualSolution{T}(blocks_from_vector(optimizer.blockdims, conv(Z)))
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

# Maps the solution of the problem written for SDPA (`data.post`) to the variables of
# the reduced problem (or of the model, in `:fallback` mode), and restores the rows
# whose constants were replaced.
function unscale_post!(optimizer::Optimizer{T}) where {T}
    data = optimizer.presolve_data
    with_work_precision(T) do
        H = eltype(data.b)
        post = data.post
        y = H.(optimizer.y)
        y[data.keep] .*= post.col ./ post.rhs
        optimizer.y = T[_to(T, data.pT, v) for v in y]
        scale_cone_values!(optimizer, post.row ./ post.obj, one(H) ./ (post.row .* post.rhs), data.pT;
                           rows = data.const_rows, values = data.const_vals)
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

# The number of `Nonnegatives` rows of `sets` (they come first, see `constraint_matrix`).
_nonnegative_rows(sets) = sum((MOI.dimension(S) for S in sets if S isa MOI.Nonnegatives); init = 0)

# The fallback problem: the cone constraints `A x + b ∈ sets` and the equality
# constraints `AE x + bE = 0` as one conic problem, with `AE x + bE ≥ 0` and
# `-AE x - bE ≥ 0` as an extra `Nonnegatives` block after the `Nonnegatives` blocks.
# Also returns the row factors `row`, `eqrow` in the same order.
function fallback_problem(sets, A, b, AE, bE, row, eqrow)
    nn = _nonnegative_rows(sets)
    psd = (nn+1):length(b)
    setsf = Any[S for S in sets if S isa MOI.Nonnegatives]
    isempty(bE) || push!(setsf, MOI.Nonnegatives(2 * length(bE)))
    append!(setsf, [S for S in sets if !(S isa MOI.Nonnegatives)])
    return setsf, vcat(A[1:nn, :], AE, -AE, A[psd, :]), vcat(b[1:nn], bE, -bE, b[psd]),
           vcat(row[1:nn], eqrow, eqrow, row[psd])
end

# Switches to the fallback and returns its model: the model without reduction, with
# each equality constraint written as two `Nonnegatives` rows (one extra block after
# the `Nonnegatives` blocks), scaled like the reduced problem.
function fallback_model!(optimizer::Optimizer{T}, reason::String) where {T}
    data = optimizer.presolve_data
    data.mode = :fallback
    data.fallback_reason = reason
    if optimizer.verbosity != SILENT
        @warn "SDPAFamily presolve: $reason. The problem is passed to the solver without the reductions of the presolve."
    end
    # for the tests, which count the fallbacks
    @debug "SDPAFamily presolve fallback" reason _group = :presolve_fallback
    # the fallback model has the original objective
    optimizer.objconstant -= data.objconstant_delta
    data.objconstant_delta = zero(T)
    data.elim = nothing
    data.redundant = Int[]
    A, b, c, AE, bE = data.original
    m = length(bE)
    nnn = count(S -> S isa MOI.Nonnegatives, data.sets)
    data.fallback_block = m == 0 ? 0 : nnn + 1
    return with_work_precision(T) do
        H = eltype(b)
        setsf, Af, bf, _ = fallback_problem(data.sets, A, b, AE, bE, ones(H, length(b)), ones(H, m))
        # replace only the constants that hold exactly; SDPA decides about the others
        crows, replacement = constant_rows(setsf, Af)
        holds = trues(length(crows))
        offset = 0
        for S in setsf
            rows = (offset+1):(offset+MOI.dimension(S))
            inblock = [i in rows for i in crows]
            if S isa MOI.Nonnegatives
                holds[inblock] .= bf[crows[inblock]] .>= 0
            elseif any(inblock)
                # a PSD block without coefficients is replaced as a whole
                holds[inblock] .= is_psd(blocks_from_vector([MOI.side_dimension(S)], bf[rows])[1])
            end
            offset = last(rows)
        end
        crows, replacement = crows[holds], replacement[holds]
        model, post, scaled = scaled_sdpa_model(T, data.pT, Af, bf, c, setsf, crows, replacement)
        data.post = post
        data.sdpa_problem = (scaled..., setsf)
        data.keep = collect(1:data.nvars)
        data.const_rows = crows
        data.const_vals = bf[crows]
        model
    end
end

# Maps the solution of the problem solved by SDPA (in `optimizer.y`, with zeros for the
# redundant variables already inserted) to the model and checks it. In `:reduced`
# mode, returns `nothing`, or the failed check (then the caller falls back). In
# `:fallback` mode, a failed check sets the status `:unverified`; returns `nothing`.
function postsolve!(optimizer::Optimizer{T}) where {T}
    data = optimizer.presolve_data
    data === nothing && return nothing
    primal_ray = MOI.get(optimizer, MOI.PrimalStatus()) == MOI.INFEASIBILITY_CERTIFICATE
    dual_ray = MOI.get(optimizer, MOI.DualStatus()) == MOI.INFEASIBILITY_CERTIFICATE
    cert_reason = data.post === nothing ? nothing : check_sdpa_certificate(optimizer)
    data.post === nothing || unscale_post!(optimizer)
    reason = cert_reason !== nothing ? cert_reason : with_work_precision(T) do
        H = eltype(data.b)
        conv(v) = T[_to(T, data.pT, vi) for vi in v]
        m = length(data.bE)
        gtol = check_tol(optimizer, data)
        unit = data.post === nothing ? one(H) : one(H) / (data.post.obj * data.post.rhs)
        if data.mode == :fallback
            _, _, _, AE, bE = data.original
            x = H.(optimizer.y)
            data.eq_primal = conv(primal_ray ? AE * x : AE * x .+ bE)
            y = weighted_dual_vector(optimizer, optimizer.X, H)
            nn = _nonnegative_rows(data.sets)
            data.eq_dual = conv(y[(nn+1):(nn+m)] .- y[(nn+m+1):(nn+2m)])
            # check the problem solved by SDPA, with the rows scaled by `data.pre`
            setsf, Af, bf, rowf = fallback_problem(data.sets, data.A, data.b, data.AE, data.bE,
                                                   data.pre.row, data.pre.eqrow)
            return check_result(optimizer, setsf, Af, bf, data.c, spzeros(H, 0, length(x)), H[],
                                x, y ./ rowf, H[], gtol, unit)
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
        all(isfinite, r) || return "the point is not finite"
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
            abs(r[i]) <= bound || return "the original equality constraints do not hold at the point"
        end
        optimizer.y = xT
        data.eq_primal = conv(r)
        yK = weighted_dual_vector(optimizer, optimizer.X, H)
        g = dual_ray ? -(transpose(data.A) * yK) : data.c - transpose(data.A) * yK
        yE = data.elim === nothing ? H[] : equality_duals(data.elim, m, g)
        if data.elim !== nothing && MOI.get(optimizer, MOI.DualStatus()) == MOI.FEASIBLE_POINT
            # The relative residual of `A' y = c` on the original data must not exceed
            # that of the reduced problem solved by SDPA by much: the mapping back is
            # exact up to rounding, and the dropped columns are consistent up to the
            # dependence tolerance. This catches ill-conditioned eliminations that
            # amplify the solver's error.
            relres(res, scale) = maximum(abs, res; init = zero(H)) /
                                 max(maximum(scale; init = zero(H)), floatmin(Float64))
            res = g .- transpose(data.AE) * yE
            scale = abs.(data.c) .+ transpose(abs.(data.A)) * abs.(yK) .+
                    transpose(abs.(data.AE)) * abs.(yE)
            res_red = data.cred .- transpose(data.Ared) * yK
            scale_red = abs.(data.cred) .+ transpose(abs.(data.Ared)) * abs.(yK)
            relres(res, scale) <= 100 * (relres(res_red, scale_red) + data.data_tol) ||
                return "the dual residual is larger than that of the reduced problem"
        end
        data.eq_dual = conv(yE)
        return check_result(optimizer, data.sets, data.A, data.b, data.c, data.AE, data.bE,
                            xr, yK, yE, gtol, unit)
    end
    if data.mode == :fallback
        if reason !== nothing
            data.sdpa_phase = optimizer.phasevalue
            data.check_failure = reason
            optimizer.phasevalue = :unverified
            if optimizer.verbosity != SILENT
                @warn "SDPAFamily: the result of SDPA did not pass the check: $reason. The status is `NUMERICAL_ERROR`."
            end
        end
        return nothing
    end
    reason === nothing && unscale_pre!(optimizer)
    return reason
end
