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
# The problem data, given in the number type `T` of the optimizer, are taken as exact
# numbers. They are converted to `BigFloat` with twice the precision of `T` plus 64 bits
# (the working precision; `Rational` data are used as they are), and the presolve
# computes in that type. There are two tolerances (see `Tolerances`):
#
# * Rounding residue. A computed number is set to zero if its terms cancelled so far
#   that its leading bit lies more than `(precision(T) + working precision) / 2` bits
#   below the leading bit of its largest term. An exact cancellation leaves a residue
#   near the working precision, far below this threshold; numbers that are exact
#   consequences of the data (such as `L - L + s` with small `s`, or a small
#   coefficient next to a large one) are far above it unless the dynamic range of
#   one computation exceeds the threshold.
# * Dependence, with `δ = 64 eps(T)` (or the `presolve_tol` option if it is larger).
#   A row of `A_E` (or a column of `A_K`) is dependent on the rows (columns) before it
#   if a combination `w` of them satisfies `|Σ_j w_j a_j| ≤ δ Σ_j |w_j| |a_j|` in every
#   entry, computed directly from the data vectors `a_j`. This is the componentwise
#   relative backward error (Oettli and Prager): relative perturbations of the data of
#   at most `δ` make the vectors exactly dependent. The criterion does not change when
#   rows or columns are scaled. The row (column) is consistent if the same
#   combination of the right-hand sides (objective coefficients) is zero in the same
#   sense; otherwise the combination is a certificate of infeasibility (a ray).
#
# The eliminations are sparse LU factorizations with threshold pivoting (see
# `eliminate`). Among the pivot candidates the one with the fewest nonzeros is taken;
# the counts are those of the data and are not updated during the elimination (a
# static tie-break, not Markowitz pivoting). After the back substitution, the reduced
# system is checked against the original rows of `A_E` at one point. If the working
# precision turns out to be insufficient (the two ways of computing a vector disagree,
# or the check fails), the presolve falls back.
#
# ## Scaling
#
# The reduced problem is scaled before it is written for SDPA (see
# `scaled_sdpa_model`). All factors are powers of two. The steps are those of conic
# solvers such as SCS, Ruiz equilibration of the rows and columns of the matrix and
# normalization of the right-hand side and of the objective, and one more between
# them: a variable with bounds gets a unit within the range that the bounds allow
# for its magnitude. The bounds are the rows with one variable, and those that follow
# from a row with several variables and the bounds of the others (the bound
# propagation of linear programming presolve). Equilibration alone leaves the unit of such
# a variable undetermined over a wide range, and its bounds can fall below SDPA's
# tolerance.
#
# SDPA's tolerances are absolute, so they apply to the scaled problem: its
# feasibility tolerance is relative to the equilibrated rows and columns, and its
# gap tolerance is relative to `max(|objective|, u)`, where `u` is the value in the
# model of an objective value of 1 in the scaled problem. Rows without coefficients
# are constants; they are checked and replaced by strictly feasible constants.
#
# The objective takes no part in the units of the variables. If the optimal value is
# small compared to the terms of the scaled objective, SDPA's gap tolerance does not
# resolve it. For that case there is a second scaling, which equilibrates the matrix
# bordered by the right-hand side and the objective. It is used if the result with
# the first scaling is rejected by a check (see `rescaled_model!`), as linear
# programming solvers solve again when the unscaled solution violates a tolerance.
#
# ## Validation and fallback
#
# * Every result that the presolve returns without calling SDPA is validated against
#   the data at `δ`: certificates of infeasibility (`A'y = 0`, `b'y < 0`), rays
#   (`A_E d = 0`, `A_K d ∈ K`, `c'd < 0`), and the point determined by the equality
#   constraints.
# * The result of SDPA is checked on the scaled problem that SDPA was given, with the
#   stopping criteria of SDPA at the square root of its tolerances (see
#   `check_sdpa_result`); certificates of infeasibility are checked for the bound
#   they prove.
# * The solution in the units of the model is checked row by row and variable by
#   variable with the same tolerance, relative to the terms of each row (see
#   `model_violations`). A row that the scaling made smaller than SDPA's tolerance
#   is not enforced by SDPA; this check finds it.
# * After the solution is mapped back, what the presolve removed is checked at the
#   returned solution: the equality constraints (in particular the dependent rows),
#   and the dual equations of the eliminated and of the dependent variables. Each
#   must hold up to the rounding to `T` and the tolerance `δ` of the decision.
#
# A result of SDPA that fails these checks is solved once more with the second
# scaling. The status is `NUMERICAL_ERROR` if that result fails too, or if it
# contradicts the first one (see `check_rescaled_result!`).
#
# The presolve falls back if the reduction itself is in doubt: the working precision
# was not sufficient, a result of the presolve fails its validation, or what the
# presolve removed does not hold at the solution. The fallback passes the problem to
# SDPA without reduction, with each equality constraint written as two inequalities
# and the same scaling.
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
    # the dependence and consistency tolerance used (see `Tolerances`)
    data_tol::H
    # cone constraints `A x + b ∈ sets` and objective `c'x` of the model, see
    # `constraint_matrix`
    A::SparseMatrixCSC{H,Int}
    b::Vector{H}
    c::Vector{H}
    # equality constraints `AE x + bE = 0`; constraint `k` has rows `eq_rows[k]`
    AE::SparseMatrixCSC{H,Int}
    bE::Vector{H}
    eq_rows::Vector{UnitRange{Int}}
    elim::Union{Nothing,EqualityElimination{H}}
    # variables of the model after the equality elimination that are fixed to zero,
    # and for each the coefficients `w` of the combination of columns that is zero
    redundant::Vector{Int}
    dependence::Dict{Int,Dict{Int,H}}
    # the scaling of the problem passed to SDPA, and its columns among the variables of
    # the reduced problem (all variables in `:fallback` mode)
    post::Union{Nothing,Scaling{H}}
    keep::Vector{Int}
    # the rows without coefficients whose constants were replaced (see
    # `scaled_sdpa_model`), and their values before the scaling
    const_rows::Vector{Int}
    const_vals::Vector{H}
    # the scaled data `(A, b, c)` and the cone sets of the problem written for SDPA
    sdpa_problem::Any
    # the same problem before the scaling, `(A, b, c, sets, crows, replacement)`, and
    # whether it was scaled with `border = true` (see `scaled_sdpa_model`)
    unscaled_problem::Any
    border::Bool
    # whether the primal point and the dual point of the last result satisfy the model
    # (see `check_model_result`), and the same for the result with the first scaling
    # if SDPA reported them as feasible and the problem was solved again
    model_feasible::Tuple{Bool,Bool}
    first_feasible::Tuple{Bool,Bool}
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

Apply the presolve described at the top of `src/presolve.jl`. Return the model to
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
        tol = Tolerances(T, pT; data_tol = optimizer.presolve_tol)
        ctol = tol.data
        n = size(A_T, 2)
        data = PresolveData{T,H}(n, pT, sets, ctol, H.(A_T), H.(b_T), H.(c_T), H.(AE_T), H.(bE_T),
                                 eq_rows, nothing, Int[], Dict{Int,Dict{Int,H}}(), nothing, Int[],
                                 Int[], H[], nothing, nothing, false, (false, false), (false, false), :reduced, 0, "", zero(T), T[], T[],
                                 :not_called, "")
        optimizer.presolve_data = data
        fallback(reason) = fallback_model!(optimizer, reason)
        A, b, c = data.A, data.b, data.c
        cmag = abs.(c)
        delta = zero(H)
        if !isempty(data.bE)
            # the number of nonzeros of each variable, to choose pivots with little fill-in
            pivot_cost = [length(nzrange(A, j)) + length(nzrange(data.AE, j)) for j in 1:n]
            # The pivot variables are chosen in the equilibrated model: a variable is
            # not preferred or avoided as a pivot because of its unit. The ranges that
            # the bounds give are applied as in `scaled_sdpa_model`, so that a variable
            # with a small range is not computed from larger ones by cancellation.
            groups, ng = _row_groups(sets, length(data.bE))
            _, colexp = equilibrate(vcat(data.AE, A), vcat(collect(1:length(data.bE)), groups), ng)
            colexp .= clamp.(colexp, magnitude_ranges(sets, A, b)...)
            elim = _eliminate_equalities(data.AE, data.bE, tol, pivot_cost, colexp)
            elim.failed && return fallback("the equality elimination did not pass its check")
            if elim.infeasible !== nothing
                y = elim.infeasible
                valid_certificate(data, zeros(H, length(b)), y, ctol) ||
                    return fallback("the certificate for inconsistent equality constraints did not pass its check")
                set_presolve_result!(optimizer, :presolve_pINF; eq_dual = y)
                return nothing
            end
            data.elim = elim
            A, b, c, cmag, delta = substitute(A, b, c, elim, tol)
        end
        # Columns that depend on the columns before them. The consistency of the
        # objective is decided relative to the terms `cmag` of the objective of the model.
        # The pivot rows are chosen in the equilibrated matrix, and among the candidates
        # the row with the fewest nonzeros (a variable with a row of its own, such as a
        # bound or an entry of a matrix variable, then needs no elimination).
        groups, ng = _row_groups(sets, 0)
        groupexp, _ = equilibrate(A, groups, ng)
        rowcount = zeros(Int, size(A, 1))
        for r in rowvals(A)
            rowcount[r] += 1
        end
        E = eliminate(A, c, tol; cmag = cmag, pivot_scale = groupexp[groups], pivot_cost = rowcount)
        E.failed && return fallback("the working precision was not sufficient to find the dependent columns")
        if E.inconsistent != 0
            ray = zeros(H, size(A, 2))
            for (j, w) in E.certificate
                ray[j] = w
            end
            d = data.elim === nothing ? ray : expand_solution(data.elim, n, ray; homogeneous = true)
            valid_ray(data, d, ctol) || return fallback("the ray of the presolve did not pass its check")
            set_presolve_result!(optimizer, :presolve_dINF; x = d)
            optimizer.primalobj = _to(T, pT, sum(data.c .* d; init = zero(H)))
            return nothing
        end
        data.redundant = E.redundant
        data.dependence = E.dependence
        keep = setdiff(1:size(A, 2), E.redundant)
        data.objconstant_delta = optimizer.objsign * _to(T, pT, delta)
        optimizer.objconstant += data.objconstant_delta
        if isempty(keep)
            return determined_point!(optimizer, data)
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
        data.keep = keep
        data.unscaled_problem = (Ak, b, c[keep], sets, crows, replacement)
        return scaled_model!(optimizer, false)
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
function determined_point!(optimizer::Optimizer{T}, data::PresolveData{T,H}) where {T,H}
    x0, s, m = equality_point(data)
    rE, mE = _apply(data.AE, x0)
    _all_zero(rE .+ data.bE, mE .+ abs.(data.bE), data.data_tol) ||
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

# Maps the solution of the problem written for SDPA (`data.post`) to the variables of
# the reduced problem (or of the model, in `:fallback` mode), and restores the rows
# whose constants were replaced.
function unscale_post!(optimizer::Optimizer{T}) where {T}
    data = optimizer.presolve_data
    with_work_precision(T) do
        H = eltype(data.b)
        post = data.post
        y = H.(optimizer.y)
        y[data.keep] .*= post.col
        optimizer.y = T[_to(T, data.pT, v) for v in y]
        scale_cone_values!(optimizer, post.row ./ post.obj, one(H) ./ post.row, data.pT;
                           rows = data.const_rows, values = data.const_vals)
        optimizer.primalobj = _to(T, data.pT, H(optimizer.primalobj) / post.obj)
        optimizer.dualobj = _to(T, data.pT, H(optimizer.dualobj) / post.obj)
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
    return optimizer
end

# The number of `Nonnegatives` rows of `sets` (they come first, see `constraint_matrix`).
_nonnegative_rows(sets) = sum((MOI.dimension(S) for S in sets if S isa MOI.Nonnegatives); init = 0)

# The fallback problem: the cone constraints `A x + b ∈ sets` and the equality
# constraints `AE x + bE = 0` as one conic problem, with `AE x + bE ≥ 0` and
# `-AE x - bE ≥ 0` as an extra `Nonnegatives` block after the `Nonnegatives` blocks.
function fallback_problem(sets, A, b, AE, bE)
    nn = _nonnegative_rows(sets)
    psd = (nn+1):length(b)
    setsf = Any[S for S in sets if S isa MOI.Nonnegatives]
    isempty(bE) || push!(setsf, MOI.Nonnegatives(2 * length(bE)))
    append!(setsf, [S for S in sets if !(S isa MOI.Nonnegatives)])
    return setsf, vcat(A[1:nn, :], AE, -AE, A[psd, :]), vcat(b[1:nn], bE, -bE, b[psd])
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
    empty!(data.dependence)
    m = length(data.bE)
    nnn = count(S -> S isa MOI.Nonnegatives, data.sets)
    data.fallback_block = m == 0 ? 0 : nnn + 1
    return with_work_precision(T) do
        setsf, Af, bf = fallback_problem(data.sets, data.A, data.b, data.AE, data.bE)
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
        data.keep = collect(1:data.nvars)
        data.const_rows = crows
        data.const_vals = bf[crows]
        data.unscaled_problem = (Af, bf, data.c, setsf, crows, replacement)
        scaled_model!(optimizer, false)
    end
end

# Scales the problem to pass to SDPA (`data.unscaled_problem`) and returns its model.
# Call inside `with_work_precision`.
function scaled_model!(optimizer::Optimizer{T}, border::Bool) where {T}
    data = optimizer.presolve_data
    A, b, c, sets, crows, replacement = data.unscaled_problem
    model, post, scaled = scaled_sdpa_model(T, data.pT, A, b, c, sets, crows, replacement; border = border)
    data.post = post
    data.sdpa_problem = (scaled..., sets)
    data.border = border
    return model
end

# If the result of SDPA was rejected by a check (status `:unverified`) and the problem
# was scaled with the first scaling, returns the model with the second scaling, to be
# solved again; otherwise `nothing`. The scaling is a heuristic and the checks decide
# whether it served; see `scaled_sdpa_model`.
function rescaled_model!(optimizer::Optimizer{T}) where {T}
    data = optimizer.presolve_data
    (optimizer.phasevalue == :unverified && !data.border && data.check_failure != OVERFLOW) || return nothing
    if optimizer.verbosity == VERBOSE
        @info "SDPAFamily presolve: the result did not pass the check ($(data.check_failure)); solving again with the second scaling."
    end
    data.first_feasible = (data.model_feasible[1] && data.sdpa_phase in (:pdOPT, :pdFEAS, :pFEAS),
                           data.model_feasible[2] && data.sdpa_phase in (:pdOPT, :pdFEAS, :dFEAS))
    data.model_feasible = (false, false)
    data.sdpa_phase = :not_called
    data.check_failure = ""
    return with_work_precision(() -> scaled_model!(optimizer, true), T)
end

# After the solve with the second scaling: a verdict of SDPA that contradicts the result
# with the first scaling is not accepted. That is the verdict "infeasible" if the first
# result had a primal point that SDPA reported as feasible and that satisfies the model
# (see `model_violations`), and the verdict "unbounded" if it had such a dual point.
function check_rescaled_result!(optimizer::Optimizer)
    data = optimizer.presolve_data
    status = MOI.get(optimizer, MOI.TerminationStatus())
    primal, dual = data.first_feasible
    if (primal && status in (MOI.INFEASIBLE, MOI.INFEASIBLE_OR_UNBOUNDED)) ||
       (dual && status in (MOI.DUAL_INFEASIBLE, MOI.INFEASIBLE_OR_UNBOUNDED))
        _set_unverified!(optimizer, "with the first scaling of the problem SDPA returned a feasible point of the model, with the second it reports `$(optimizer.phasevalue)`")
    end
    return optimizer
end

# Checks the result of SDPA on the scaled problem, maps it back to the model (it is in
# `optimizer.y`, with zeros for the redundant variables already inserted), and checks
# what the presolve removed. A result of SDPA that fails the check on the scaled
# problem, or overflows in `T`, gets the status `:unverified`. In `:reduced` mode,
# returns the failed check of what the presolve removed (then the caller falls back),
# or `nothing`.
function postsolve!(optimizer::Optimizer{T}) where {T}
    data = optimizer.presolve_data
    (data === nothing || data.post === nothing) && return nothing
    primal_ray = MOI.get(optimizer, MOI.PrimalStatus()) == MOI.INFEASIBILITY_CERTIFICATE
    dual_ray = MOI.get(optimizer, MOI.DualStatus()) == MOI.INFEASIBILITY_CERTIFICATE
    has_dual = MOI.get(optimizer, MOI.DualStatus()) in (MOI.FEASIBLE_POINT, MOI.INFEASIBILITY_CERTIFICATE)
    reason = check_sdpa_result(optimizer)
    if optimizer.phasevalue == :pdFEAS && (reason == GAP_NOT_SMALL || reason == GAP_NOT_RESOLVED)
        # feasible iterates at the iteration limit whose gap is not small: not optimal
        optimizer.phasevalue = :pdFEAS_gap
        reason = nothing
    end
    unscale_post!(optimizer)
    if data.mode == :fallback
        with_work_precision(T) do
            H = eltype(data.b)
            conv(v) = T[_to(T, data.pT, vi) for vi in v]
            m = length(data.bE)
            x = H.(optimizer.y)
            data.eq_primal = conv(primal_ray ? data.AE * x : data.AE * x .+ data.bE)
            y = weighted_dual_vector(optimizer, optimizer.X, H)
            nn = _nonnegative_rows(data.sets)
            data.eq_dual = conv(y[(nn+1):(nn+m)] .- y[(nn+m+1):(nn+2m)])
        end
        reason === nothing && !_finite_result(optimizer) && (reason = OVERFLOW)
        model_reason = with_work_precision(T) do
            H = eltype(data.b)
            nn = _nonnegative_rows(data.sets)
            m = length(data.bE)
            y = weighted_dual_vector(optimizer, optimizer.X, H)
            yK = vcat(y[1:nn], y[(nn+2m+1):end])
            check_model_result(optimizer, H.(optimizer.y), yK, H.(data.eq_dual))
        end
        reason === nothing && (reason = model_reason)
        reason === nothing || _set_unverified!(optimizer, reason)
        return nothing
    end
    return with_work_precision(T) do
        H = eltype(data.b)
        conv(v) = T[_to(T, data.pT, vi) for vi in v]
        elim = data.elim
        # the values of all variables and of the equality constraints, and the duals of
        # the equality constraints
        x = elim === nothing ? H.(optimizer.y) :
            expand_solution(elim, data.nvars, H.(optimizer.y); homogeneous = primal_ray)
        optimizer.y = conv(x)
        xr = H.(optimizer.y)
        r, rm = _apply(data.AE, xr)
        if !primal_ray
            r .+= data.bE
            rm .+= abs.(data.bE)
        end
        data.eq_primal = conv(r)
        yK = weighted_dual_vector(optimizer, optimizer.X, H)
        g = dual_ray ? -(transpose(data.A) * yK) : data.c - transpose(data.A) * yK
        gm = transpose(abs.(data.A)) * abs.(yK)
        dual_ray || (gm .+= abs.(data.c))
        data.eq_dual = conv(elim === nothing ? H[] : equality_duals(elim, length(data.bE), g))
        # The result of SDPA for the reduced problem failed its check, or it passed and
        # overflows in `T`. The reduced problem is equivalent to the model and better
        # posed than the model with its equality constraints as pairs of
        # inequalities, so the fallback would not give a more reliable result.
        if reason === nothing && (!_finite_result(optimizer) || !all(isfinite, data.eq_dual))
            reason = OVERFLOW
        end
        model_reason = check_model_result(optimizer, xr, yK, H.(data.eq_dual))
        reason === nothing && (reason = model_reason)
        if reason !== nothing
            _set_unverified!(optimizer, reason)
            return nothing
        end
        # What the presolve removed is checked at the returned solution, relative to
        # the terms of each row (column), up to `ptol` for the rounding to `T`.
        ptol = 64 * eps_of(T, data.pT)
        δ = data.data_tol
        # 1. The equality constraints at the returned point. The rows that were kept
        # hold up to the rounding of `x` to `T`. A row `i` dropped as dependent holds up
        # to the dependence tolerance times the terms of its combination `w` with other
        # rows, plus the residuals of those rows.
        _removed_hold(r, rm, elim === nothing ? Dict{Int,Dict{Int,H}}() : elim.dependence, δ, ptol) ||
            return "the original equality constraints do not hold at the point"
        # 2. The dual equations `c - A'y_K - A_E'y_E = 0` (for a certificate, without
        # `c`) of the variables that the presolve removed. Those of the pivot variables
        # of the equality elimination hold by the construction of `y_E`. That of a
        # variable `j` dropped as dependent holds up to the tolerance of the decision
        # and the residuals of the variables in its combination.
        if has_dual
            yE = H.(data.eq_dual)
            res = g .- transpose(data.AE) * yE
            mag = gm .+ transpose(abs.(data.AE)) * abs.(yE)
            free = elim === nothing ? collect(1:data.nvars) : elim.free
            if elim !== nothing
                all(abs(res[p]) <= ptol * mag[p] for p in elim.pivots) ||
                    return "the dual equations of the eliminated variables do not hold"
                # the terms of the dual equation of a free variable in the reduced
                # problem include those of the pivot variables it replaces
                mag = copy(mag)
                for (q, p) in enumerate(elim.pivots), (j, w) in _pairs(elim.rows[q])
                    mag[j] += abs(w) * mag[p]
                end
            end
            dependence = Dict{Int,Dict{Int,H}}(free[k] => Dict{Int,H}(free[l] => wl for (l, wl) in w)
                                               for (k, w) in data.dependence)
            _removed_hold(res, mag, dependence, δ, ptol; only = keys(dependence)) ||
                return "the dual equations of the variables dropped as dependent do not hold"
        end
        return nothing
    end
end

const OVERFLOW = "the solution is not finite in the number type of the optimizer"

# Sets the status `:unverified` (`NUMERICAL_ERROR`): the result of SDPA is kept, with
# the status reported by SDPA and the failed check in `MOI.RawStatusString`.
function _set_unverified!(optimizer::Optimizer, reason::String)
    data = optimizer.presolve_data
    data.sdpa_phase = optimizer.phasevalue
    data.check_failure = reason
    optimizer.phasevalue = :unverified
    if optimizer.verbosity != SILENT
        @warn "SDPAFamily: the result of SDPA did not pass the check: $reason. The status is `NUMERICAL_ERROR`."
    end
    return optimizer
end

# Whether the solution stored in `optimizer` (in the number type `T`) is finite. The
# scaled solution of SDPA is; the solution of the model can overflow in `T`.
_finite_result(optimizer::Optimizer) =
    all(isfinite, optimizer.y) && isfinite(optimizer.primalobj) && isfinite(optimizer.dualobj) &&
    all(B -> all(isfinite, B), optimizer.X.blocks) && all(B -> all(isfinite, B), optimizer.Z.blocks)

# Whether the residuals `r` (with the sums `rm` of the absolute values of their terms)
# of the rows `only` are zero up to `ptol`, or, for a row `i` with a dependence `w`
# (`w[i] = 1`, see `Elimination`), up to the tolerance `δ` of the dependence decision
# and the residuals of the other rows of `w`.
function _removed_hold(r::Vector{H}, rm::Vector{H}, dependence, δ, ptol; only = eachindex(r)) where {H}
    for i in only
        bound = ptol * rm[i]
        if haskey(dependence, i)
            for (k, w) in dependence[i]
                k == i && continue
                bound += abs(w) * (2 * δ * rm[k] + abs(r[k]))
            end
            bound += 2 * δ * rm[i]
        end
        abs(r[i]) <= bound || return false
    end
    return true
end
