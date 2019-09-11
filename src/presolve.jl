using SparseArrays
export presolve

"""
    presolve(optimizer::SDPAFamily.Optimizer{T}) where T

Identifies linearly dependent constraints in the problem. This is done by a naive Gaussian elimination.

Returns a vector with the indices of redundant constraints, which should be removed from the formulation. The rest of the constraints form a maximal linearly independent subset of the original set of constraints.
"""
function presolve(optimizer::SDPAFamily.Optimizer{T}) where T
    start = time()
    totaldim = sum(abs.(optimizer.blockdims))
    F = spzeros(T, length(optimizer.b), totaldim^2)
    cVec = optimizer.b
    for entry in optimizer.elemdata
        constr_index, blk, i, j, value = entry
        if constr_index != 0
            offset = sum(abs.(optimizer.blockdims)[1:blk-1])
            linear_index = (j + offset - 1) * totaldim + (i + offset)
            F[constr_index, linear_index] = value
        end
    end
    aug_mat = hcat(F, cVec)
    redundant_F = collect(setdiff!(
        Set(1:length(optimizer.b)),
        Set(rowvals(reduce!(aug_mat)[:, 1:end-1]))
    ))
    abort = 0
    for i in redundant_F
        if aug_mat[i, end] != 0
            # println(aug_mat[i, end])
            abort = 1
            @warn "Inconsistency at constraint index $i. Problem is dual infeasible."
        end
    end
    if abort == 1
        optimizer.phasevalue = :pFEAS_dINF
        optimizer.y = zeros(T, length(optimizer.b))
        optimizer.X = PrimalSolution{T}(map(
            n -> zeros(T, abs(n), abs(n)),
            optimizer.blockdims
        ))
        optimizer.Z = VarDualSolution{T}(map(
            n -> zeros(T, abs(n), abs(n)),
            optimizer.blockdims
        ))
    end
    finish = time() - start
    n = length(redundant_F)
    if !optimizer.silent
        @info "Presolve finished in $finish seconds. $n constraint(s) eliminated."
    end
    return sort!(redundant_F)
end

"""
    function reduce!(A::SparseMatrixCSC{T}, ɛ = T <: Union{Rational,Integer} ? 0 : eps(norm(A, Inf))) where T

Identifies linearly dependent constraints in the problem. The last column of input is constraint constants and they are included to check if the linearly dependent constraints are redundant or inconsistent. This is done by a naive Gaussian elimination.

Returns a vector with the indices of redundant constraints, which should be removed from the formulation. The rest of the constraints form a maximal linearly independent subset of the original set of constraints.
"""
function reduce!(A::SparseMatrixCSC{T}, ɛ = T <: Union{Rational,Integer} ? 0 : eps(norm(A, Inf))) where T
    nr, nc = size(A)
    i = j = 1
    visited_rows = Int[]
    while i <= nr && j <= nc - 1 # avoid touching the cVec
        rows = setdiff!(rowvals(A)[nzrange(A, j)], visited_rows)
        if isempty(rows)
            j += 1
            continue
        end
        (m, mi) = findmax(abs.([A[p, j] for p in rows]))
        mi = rows[mi]
        if m <= ɛ
            if ɛ > 0
                A[:, j] .= zero(T)
                droptol!(A, ɛ)
            end
            j += 1
        else
            push!(visited_rows, mi)
            d = A[mi, j]
            A[mi, j:nc] = A[mi, j:nc] / d

            for k in setdiff!(rows, mi)
                d = A[k, j]
                A[k, j:nc] -= d * A[mi, j:nc]
            end
            droptol!(A, ɛ)
            i += 1
            j += 1
        end
    end
    return droptol!(A, ɛ)
end
