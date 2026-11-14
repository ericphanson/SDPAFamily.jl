# Random feasible and bounded conic problems (the generator of the review of the
# presolve): a strictly feasible primal point `x0` and a strictly feasible dual point
# give the constraints and the objective. The problems have redundant and duplicated
# equality constraints, duplicated variable columns, and (with `scale = true`) equality
# rows scaled by 10^±12 and `Nonnegatives` rows by 10^±6, which does not change the
# solution.
function random_conic_problem(rng, ::Type{T}; scale::Bool, nmax = 10, emax = 6, nnmax = 4,
                              dmax = 5) where {T}
    n = rand(rng, 2:nmax)
    floatdata = rand(rng) < 0.7
    gen() = floatdata ? T(round(randn(rng), digits = 1)) : T(rand(rng, -3:3))
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
    if scale
        for i in 1:size(AE, 1)
            s = T(10)^rand(rng, (-12, 0, 12))
            AE[i, :] .*= s
            bE[i] *= s
        end
        for i in 1:size(AN, 1)
            s = T(10)^rand(rng, (-6, 0, 6))
            AN[i, :] .*= s
            bN[i] *= s
        end
    end
    sense = rand(rng, (MOI.MIN_SENSE, MOI.MAX_SENSE))
    sense == MOI.MAX_SENSE && (c = -c)
    return (; AE, bE, AN, bN, AP, bP = d > 0 ? bP : T[], d, c, sense)
end

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
