export read_results!

"""
    read_results!(optimizer::SDPAGMPOptimizer{T}, filepath::String)
Populates `optimizer` with results in a SDPA-formatted output file specified by `filepath`.

"""
function read_results!(optimizer::SDPAGMPOptimizer{T}, filepath::String) where T

    endswith(filepath, ".dat") || error("Filename '$filepath' must end with .dat")

    getnextline(io::IO) = eof(io) ? error("The output file is possibly corrupted. Check that $filepath conforms to the SDPA output format.") : chomp(readline(io))
    
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

    phasevalue = "noINFO"
    objValPrimalstring = ""
    objValDualstring = ""
    xVecstring = ""
    xMatvec = T[]
    yMatvec = T[]

    open(filepath, "r") do io
        line = getnextline(io)

        while !startswith(line, "phase.value")
            line = getnextline(io)
        end
        phasevalue = split(line)[3]

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
    # println(stdout, phasevalue)
    xVecstring = remove_brackets!(xVecstring)
    xVecstring = split(xVecstring, ",")
    xVec = parse.(T, xVecstring)
    objValPrimal = parse(T, objValPrimalstring)
    objValDual = parse(T, objValDualstring)

    # xMatstring = replace(remove_brackets(xMatstring), " " => "")
    # yMatstring = replace(remove_brackets(yMatstring), " " => "")
    #
    # xMatvec = parse.(T, split(xMatstring[2:end], ","))
    # yMatvec = parse.(T, split(yMatstring[2:end], ","))
    

    if phasevalue == "noINFO"
        optimizer.terminationstatus = MOI.OPTIMIZE_NOT_CALLED
        optimizer.primalstatus = MOI.UNKNOWN_RESULT_STATUS
        optimizer.dualstatus = MOI.UNKNOWN_RESULT_STATUS
    elseif phasevalue == "pFEAS"
        optimizer.terminationstatus = MOI.SLOW_PROGRESS
        optimizer.primalstatus = MOI.FEASIBLE_POINT
        optimizer.dualstatus = MOI.UNKNOWN_RESULT_STATUS
    elseif phasevalue == "dFEAS"
        optimizer.terminationstatus = MOI.SLOW_PROGRESS
        optimizer.primalstatus = MOI.UNKNOWN_RESULT_STATUS
        optimizer.dualstatus = MOI.FEASIBLE_POINT
    elseif phasevalue == "pdFEAS"
        optimizer.terminationstatus = MOI.OPTIMAL
        optimizer.primalstatus = MOI.FEASIBLE_POINT
        optimizer.dualstatus = MOI.FEASIBLE_POINT
    elseif phasevalue == "pdINF"
        optimizer.terminationstatus = MOI.INFEASIBLE_OR_UNBOUNDED
        optimizer.primalstatus = MOI.UNKNOWN_RESULT_STATUS
        optimizer.dualstatus = MOI.UNKNOWN_RESULT_STATUS
    elseif phasevalue == "pFEAS_dINF"
        optimizer.terminationstatus = MOI.DUAL_INFEASIBLE
        optimizer.primalstatus = MOI.INFEASIBILITY_CERTIFICATE
        optimizer.dualstatus = MOI.INFEASIBLE_POINT
    elseif phasevalue == "pINF_dFEAS"
        optimizer.terminationstatus = MOI.INFEASIBLE
        optimizer.primalstatus = MOI.INFEASIBLE_POINT
        optimizer.dualstatus = MOI.INFEASIBILITY_CERTIFICATE
    elseif phasevalue == "pdOPT"
        optimizer.terminationstatus = MOI.OPTIMAL
        optimizer.primalstatus = MOI.FEASIBLE_POINT
        optimizer.dualstatus = MOI.FEASIBLE_POINT
    elseif phasevalue == "pUNBD"
        optimizer.terminationstatus = MOI.DUAL_INFEASIBLE
        optimizer.primalstatus = MOI.INFEASIBILITY_CERTIFICATE
        optimizer.dualstatus = MOI.INFEASIBLE_POINT
    elseif phasevalue == "dUNBD"
        optimizer.terminationstatus = MOI.INFEASIBLE
        optimizer.primalstatus = MOI.INFEASIBLE_POINT
        optimizer.dualstatus = MOI.INFEASIBILITY_CERTIFICATE
    end

    optimizer.y = xVec

    inputpath = replace(filepath, "output.dat" => "input.dat-s")
    structurevec = []
    open(inputpath, "r") do io
        line = getnextline(io)
        line = getnextline(io)
        line = getnextline(io)
        structurevec = parse.(Int, split(line))
    end
    xMatbm = BlockMatrix{T}(map(n -> zeros(T, abs(n), abs(n)), optimizer.blkdims))
    yMatbm = BlockMatrix{T}(map(n -> zeros(T, abs(n), abs(n)), optimizer.blkdims))
    for i in 1:length(structurevec)
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
