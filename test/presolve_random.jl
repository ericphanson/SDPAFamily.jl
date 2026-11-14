# Random feasible and bounded conic problems (the generator of the review of the
# presolve): a strictly feasible primal point `x0` and a strictly feasible dual point
# give the constraints and the objective. The problems have redundant and duplicated
# equality constraints, duplicated variable columns, and (with `scale = true`) equality
# rows scaled by 10^±12 and `Nonnegatives` rows by 10^±6, which does not change the
# solution.
#
# With `kind = :infeasible`, the constraints `a'x + β >= 0` and `-a'x - β - 1 >= 0` are
# added. With `kind = :unbounded`, the constraints are changed so that a direction `d`
# is a ray (`AE d = 0`, `AN d >= 0`, `AP d` positive definite), and the objective so
# that `c'd < 0`. For these kinds the matrices are multiples of 1/16 and the row
# factors are powers of two, so that these properties hold exactly.
function random_conic_problem(rng, ::Type{T}; scale::Bool, kind::Symbol = :optimal, nmax = 10,
                              emax = 6, nnmax = 4, dmax = 5) where {T}
    n = rand(rng, 2:nmax)
    floatdata = rand(rng) < 0.7
    gen() = !floatdata ? T(rand(rng, -3:3)) :
            kind == :optimal ? T(round(randn(rng), digits = 1)) : T(round(4 * randn(rng)) / 4)
    nE = rand(rng, 0:emax)
    AE = [gen() for _ in 1:nE, _ in 1:n]
    if nE >= 1 && rand(rng) < 0.6
        λ = [gen() for _ in 1:nE]
        AE = vcat(AE, λ' * AE)
        rand(rng) < 0.5 && (AE = vcat(AE, AE[1:1, :]))
    end
    nNN = rand(rng, 1:nnmax)
    AN = [gen() for _ in 1:nNN, _ in 1:n]
    d = rand(rng, 0:dmax)
    AP = [gen() for _ in 1:div(d * (d + 1), 2), _ in 1:n]
    if n >= 3 && rand(rng) < 0.5
        AE[:, n] = AE[:, 1]
        AN[:, n] = AN[:, 1]
        AP[:, n] = AP[:, 1]
    end
    if kind == :unbounded
        dir = T.(rand(rng, (-2, -1, 1, 2), n))
        dd = sum(abs2, dir)
        AE = (dd * AE - (AE * dir) * dir') / 4
        for i in 1:nNN
            dot(AN[i, :], dir) < 0 && (AN[i, :] .*= -1)
        end
        if d > 0
            P = [T(i == j ? 2 : rand(rng, -1:1) / 4) for i in 1:d, j in 1:d]
            P = (P + P') / 2
            AP = (dd * AP + ([P[i, j] for j in 1:d for i in 1:j] * dd - AP * dir) * dir') / 4
        end
    end
    x0 = [gen() for _ in 1:n]
    bE = -AE * x0
    bN = T.(rand(rng, 1:3, nNN)) - AN * x0
    c = AN' * T.(rand(rng, 1:3, nNN)) + AE' * [gen() for _ in 1:size(AE, 1)]
    if d > 0
        S = T.(randn(rng, d, d))
        S = S * S' + I
        bP = [S[i, j] for j in 1:d for i in 1:j] - AP * x0
        Y = T.(randn(rng, d, d))
        Y = Y * Y' + I
        c += AP' * [(i == j ? 1 : 2) * Y[i, j] for j in 1:d for i in 1:j]
    end
    if kind == :unbounded
        c = dd * c - dir * (dd + dot(c, dir))
    elseif kind == :infeasible
        a = [gen() for _ in 1:n]
        iszero(a) && (a[1] = one(T))
        β = gen()
        AN = vcat(AN, a', -a')
        bN = vcat(bN, β, -β - 1)
    end
    if scale
        for i in 1:size(AE, 1)
            s = (kind == :optimal ? T(10) : T(2))^rand(rng, (-12, 0, 12))
            AE[i, :] .*= s
            bE[i] *= s
        end
        for i in 1:size(AN, 1)
            s = (kind == :optimal ? T(10) : T(2))^rand(rng, (-6, 0, 6))
            AN[i, :] .*= s
            bN[i] *= s
        end
    end
    sense = rand(rng, (MOI.MIN_SENSE, MOI.MAX_SENSE))
    sense == MOI.MAX_SENSE && (c = -c)
    return (; AE, bE, AN, bN, AP, bP = d > 0 ? bP : T[], d, c, sense)
end

# The problem `p` with variable `j` replaced by `s[j]` times a new variable.
scale_variables(p::NamedTuple, s::Vector) =
    merge(p, (; AE = p.AE .* s', AN = p.AN .* s', AP = p.AP .* s', c = p.c .* s))

function add_random_conic_problem!(model, p::NamedTuple, ::Type{T}) where {T}
    n = length(p.c)
    x = MOI.add_variables(model, n)
    function add(A, b, S)
        terms = [MOI.VectorAffineTerm(i, MOI.ScalarAffineTerm(A[i, j], x[j]))
                 for i in 1:size(A, 1) for j in 1:n if !iszero(A[i, j])]
        MOI.add_constraint(model, MOI.VectorAffineFunction(terms, b), S)
    end
    isempty(p.bE) || add(p.AE, p.bE, MOI.Zeros(length(p.bE)))
    add(p.AN, p.bN, MOI.Nonnegatives(length(p.bN)))
    p.d > 0 && add(p.AP, p.bP, MOI.PositiveSemidefiniteConeTriangle(p.d))
    MOI.set(model, MOI.ObjectiveSense(), p.sense)
    MOI.set(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}(),
            MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(p.c, x), zero(T)))
    return x
end
