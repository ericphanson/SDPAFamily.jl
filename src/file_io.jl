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
    # FIXME `optimizer.b` has been removed
    #if optimizer.verbosity != SILENT && (norm(optimizer.b, Inf) < eps(norm(xMatvec, Inf)) || norm(optimizer.b, Inf) < eps(norm(yMatvec, Inf)))
    #    @warn "Potential underflow detected. Check the results and use `BigFloat` entries if necessary."
    #end
    xVecstring = remove_brackets!(xVecstring)
    xVecstring = split(xVecstring, ",")
    xVec = parse.(T, xVecstring)
    # The variables of the MOI model are SDPA's `xVec`, so SDPA's primal objective
    # `objValPrimal = c'x` is the primal objective of the MOI model, and
    # `objValDual = F_0 • Y` is its dual objective.
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
