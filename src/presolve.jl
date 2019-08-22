using BlockArrays
using SparseArrays
using GenericLinearAlgebra
export presolve


function presolve(optimizer::SDPA_GMP.Optimizer{T}) where T
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
    F = Array((hcat(vec.(F)...))')
    SVD = svd!(F)
    S = filter(x -> abs(x)>=1e-10, SVD.S)
    U = SVD.U[:, 1:length(S)]
    # Vt = SVD.Vt[1:length(S), :]
    cVec = U'*cVec
    # F = Diagonal(S)*SVD.Vt
    # F0 = reshape(SVD.Vt*vec(Array(F0)), totaldim, totaldim)
    open(joinpath(optimizer.tempfile, "F.dat"),"w") do io
        for i in 1:size(F,1)
        end
    end
    return size(SVD.Vt)
end




