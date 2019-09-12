export read_results!

"""
    read_results!(optimizer::Optimizer{T}, filepath::String, redundant_entries::Vector)

Populates `optimizer` with results in a SDPA-formatted output file specified by `filepath`. Redundant entries corresponding to linearly dependent constraints are set to 0.
"""
function read_results!(
    optimizer::Optimizer{T},
    filepath::String,
    redundant_entries::Vector
) where T

    endswith(filepath, ".dat") || error("Filename '$filepath' must end with .dat")
    getnextline(io::IO) =
        eof(io) ?
        error("The output file is possibly corrupted. Check that $filepath conforms to the SDPA output format.") :
        chomp(readline(io))


    function replace_brackets!(str::SubString)
        str = replace(str, "{" => "[")
        str = replace(str, "}" => "]")
        return str
    end
    function remove_brackets!(str)
        str = replace(str, "{" => "")
        str = replace(str, "}" => "")
        return str
    end

    phasevalue = :noINFO
    objValPrimalstring = ""
    objValDualstring = ""
    xVecstring = ""
    xMatvec = T[]
    yMatvec = T[]

    file = open(filepath, "r") do io
            line = getnextline(io)

            while !startswith(line, "phase.value")
                line = getnextline(io)
            end
            optimizer.phasevalue = Symbol(split(line)[3])

            while !startswith(line, "objValPrimal")
                line = getnextline(io)
            end
            objValPrimalstring = split(line)[3]

            while !startswith(line, "objValDual")
                line = getnextline(io)
            end
            objValDualstring = split(line)[3]

            while !startswith(line, "xVec")
                line = getnextline(io)
            end
            line = getnextline(io)
            xVecstring = line

            while !startswith(line, "xMat")
                line = getnextline(io)
            end
            line = getnextline(io)
            line = getnextline(io)
            while !startswith(line, "}")
                xMatstring = remove_brackets!(line)
                if endswith(xMatstring, ",")
                    append!(xMatvec, parse.(T, split(xMatstring[1:end-1], ",")))
                else
                    append!(xMatvec, parse.(T, split(xMatstring, ",")))
                end
                line = getnextline(io)
            end
            line = getnextline(io)
            line = getnextline(io)
            line = getnextline(io)
            while !startswith(line, "}")
                yMatstring = remove_brackets!(line)
                if endswith(yMatstring, ",")
                    append!(yMatvec, parse.(T, split(yMatstring[1:end-1], ",")))
                else
                    append!(yMatvec, parse.(T, split(yMatstring, ",")))
                end
                line = getnextline(io)
            end
        end
    if norm(optimizer.b, Inf) < eps(norm(xMatvec, Inf)) || norm(optimizer.b, Inf) < eps(norm(yMatvec, Inf))
        @warn "Potential underflow detected. Check the results and use `BigFloat` entries if necessary."
    end
    xVecstring = remove_brackets!(xVecstring)
    xVecstring = split(xVecstring, ",")
    xVec = parse.(T, xVecstring)
    optimizer.primalobj = parse(T, objValPrimalstring)
    optimizer.dualobj = parse(T, objValDualstring)
    for i in redundant_entries
        splice!(xVec, i:i-1, zero(T))
    end
    optimizer.y = xVec
    structurevec = optimizer.blockdims
    yMatbm = PrimalSolution{T}(map(
        n -> zeros(T, abs(n), abs(n)),
        optimizer.blockdims
    ))
    xMatbm = VarDualSolution{T}(map(
        n -> zeros(T, abs(n), abs(n)),
        optimizer.blockdims
    ))
    for i = 1:length(structurevec)
        dim = structurevec[i]
        if dim < 0
            dim = abs(dim)
            xblock = xMatvec[1:dim]
            deleteat!(xMatvec, 1:dim)
            xblock = Diagonal(xblock)
            xMatbm.blocks[i] = xblock

            yblock = yMatvec[1:dim]
            deleteat!(yMatvec, 1:dim)
            yblock = Diagonal(yblock)
            yMatbm.blocks[i] = yblock


        elseif dim > 0
            xblock = xMatvec[1:dim^2]
            deleteat!(xMatvec, 1:dim^2)
            xblock = reshape(xblock, dim, dim)
            xMatbm.blocks[i] = xblock

            yblock = yMatvec[1:dim^2]
            deleteat!(yMatvec, 1:dim^2)
            yblock = reshape(yblock, dim, dim)
            yMatbm.blocks[i] = yblock

        end
    end
    optimizer.X = yMatbm
    optimizer.Z = xMatbm

end

"""
    inputElement(optimizer::Optimizer, constr_number::Int, blk::Int, i::Int, j::Int, value::T) where T

Stores the constraint data in `optimizer.elemdata` as a vector of tuples. Each tuple corresponds to one line in the SDPA-formatted input file.
"""
function inputElement(
    optimizer::Optimizer,
    constr_number::Int,
    blk::Int,
    i::Int,
    j::Int,
    value::T
) where T
    push!(optimizer.elemdata, (constr_number, blk, i, j, value))
end

"""
    initializeSolve(optimizer::Optimizer)

Writes problem data into an SDPA-formatted file named `input.dat-s`. `presolve.jl` routine is applied as indicated by `optimizer.presolve`.

Returns a vector of indices for redundant constraints, which are omitted from the input file.
"""
function initializeSolve(optimizer::Optimizer)
    if !optimizer.presolve
        filename = joinpath(optimizer.tempdir, "input.dat-s")
        file = open(filename, "w") do io
            nconstrs = length(optimizer.b)
            nblocks = length(optimizer.blockdims)
            println(io, nconstrs)
            println(io, nblocks)
            writedlm(io, optimizer.blockdims', " ")
            writedlm(io, optimizer.b', " ")
            for line in optimizer.elemdata
                for i in line
                    print(io, i)
                    print(io, " ")
                end
                println(io)
            end
        end
        return []
    else
        redundant_F = presolve(optimizer)
        if optimizer.phasevalue != :pFEAS_dINF
            reduced = joinpath(optimizer.tempdir, "input.dat-s")
            file = open(reduced, "w") do io
                nconstrs = length(optimizer.b) - length(redundant_F)
                nblocks = length(optimizer.blockdims)
                println(io, nconstrs)
                println(io, nblocks)
                writedlm(io, optimizer.blockdims', " ")
                cVec = deleteat!(copy(optimizer.b), redundant_F)
                writedlm(io, cVec', " ")
                for entry in optimizer.elemdata
                    constr_index, blk, i, j, value = entry
                    if !in(constr_index, redundant_F)
                        constr_index = constr_index -
                                       count(x -> x < constr_index, redundant_F)
                        for i in (constr_index, blk, i, j, value)
                            print(io, i)
                            print(io, " ")
                        end
                        println(io)
                    end
                end
            end
        end
        return redundant_F
    end
end
