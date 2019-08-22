using BlockArrays
using SparseArrays
using GenericLinearAlgebra
using RowEchelon
export presolve


function presolve(optimizer::SDPA_GMP.Optimizer{T}) where T
    start = time()
    filepath = joinpath(optimizer.tempfile, "input.dat-s")
    totaldim = sum(abs.(optimizer.blockdims))
    F0 = BlockArray(spzeros(T, totaldim, totaldim), abs.(optimizer.blockdims), abs.(optimizer.blockdims))
    F = Vector{Any}()
    for i in 1:length(optimizer.b)
        push!(F,BlockArray(spzeros(T, totaldim, totaldim), abs.(optimizer.blockdims), abs.(optimizer.blockdims)))
    end
    cVec = T[]
    file = open(filepath, "r") do io
        chomp(readline(io))
        chomp(readline(io))
        chomp(readline(io))
        cVec = parse.(T, split(strip(chomp(readline(io))), " "))
        while !eof(io)
            constr_index, blk, i, j, value = parse.(T, split(strip(chomp(readline(io))), " "))
            constr_index = Int(constr_index)
            blk = Int(blk)
            i = Int(i)
            j = Int(j)
            if constr_index == 0
                F0.blocks[blk, blk][i,j] = value
            else
                F[constr_index].blocks[blk, blk][i,j] = value
            end
        end
    end
    redundant_F = []
    scratch_matrix = zeros(T, 1, totaldim^2)
    oldpivot = Int64[]
    for i in 1:length(optimizer.b)
        # show(vec(F[i])')
        scratch_matrix = [scratch_matrix; vec(F[i])']
        $, pivot = rref_with_pivots!(scratch_matrix)
        if length(oldpivot) == length(pivot)
            push!(redundant_F, i)
        elseif length(oldpivot) == length(pivot) - 1
            oldpivot = pivot
        else
            error("Reduced echelon form failed: this should never happen...")
        end
    end


    # F = Array((hcat(vec.(F)...))')
    # SVD = svd!(F)
    # S = filter(x -> abs(x)>=1e-10, SVD.S)
    # U = SVD.U[:, 1:length(S)]
    # # Vt = SVD.Vt[1:length(S), :]
    # cVec = U'*cVec
    # F = Diagonal(S)*SVD.Vt
    # F0 = reshape(SVD.Vt*vec(Array(F0)), totaldim, totaldim)
    # A, pivot = rref_with_pivots(F)
    finish = time() - start
    n = length(redundant_F)
    @info "Presolve finished in $finish seconds. $n constraints are eliminated."
    return redundant_F
end
