using MathOptInterface
using SemidefiniteOptInterface
const MOI = MathOptInterface
const SDOI = SemidefiniteOptInterface
import SemidefiniteOptInterface.block
# export read_results

nextline(io::IO) = chomp(readline(io))

function replace_brackets(str::SubString)
    str = replace(str, "{" => "[")
    str = replace(str, "}" => "]")
    return str
end
function remove_brackets(str)
    str = replace(str, "{" => "")
    str = replace(str, "}" => "")
    return str
end

function read_results(optimizer, filepath::String)
    endswith(filepath, ".dat") || error("Filename '$filepath' must end with .dat")

    phasevalue = "noINFO"
    objValPrimalstring = ""
    objValDualstring = ""
    xVecstring = ""
    xMatvec = BigFloat[]
    yMatvec = BigFloat[]

    open(filepath, "r") do io
        line = nextline(io)

        while !startswith(line, "phase.value")
            line = nextline(io)
        end
        phasevalue = split(line)[3]

        while !startswith(line, "objValPrimal")
            line = nextline(io)
        end
        objValPrimalstring = split(line)[3]

        while !startswith(line, "objValDual")
            line = nextline(io)
        end
        objValDualstring = split(line)[3]

        while !startswith(line, "xVec")
            line = nextline(io)
        end
        line = nextline(io)
        xVecstring = line

        while !startswith(line, "xMat")
            line = nextline(io)
        end
        line = nextline(io)
        line = nextline(io)
        while !startswith(line, "}")
            xMatstring = remove_brackets(line)
            if endswith(xMatstring, ",")
                show(xMatstring[1:end-1])
                append!(xMatvec, parse.(BigFloat, split(xMatstring[1:end-1], ",")))
            else
                append!(xMatvec, parse.(BigFloat, split(xMatstring, ",")))
            end
            line = nextline(io)
        end
        line = nextline(io)
        line = nextline(io)
        line = nextline(io)
        while !startswith(line, "}")
            yMatstring = remove_brackets(line)
            if endswith(yMatstring, ",")
                show(yMatstring[1:end-1])
                append!(yMatvec, parse.(BigFloat, split(yMatstring[1:end-1], ",")))
            else
                append!(yMatvec, parse.(BigFloat, split(yMatstring, ",")))
            end
            line = nextline(io)
        end
    end
    # println(stdout, phasevalue)
    xVecstring = remove_brackets(xVecstring)
    xVecstring = split(xVecstring, ",")
    xVec = parse.(BigFloat, xVecstring)
    objValPrimal = parse(BigFloat, objValPrimalstring)
    objValDual = parse(BigFloat, objValDualstring)

    # xMatstring = replace(remove_brackets(xMatstring), " " => "")
    # yMatstring = replace(remove_brackets(yMatstring), " " => "")
    #
    # xMatvec = parse.(BigFloat, split(xMatstring[2:end], ","))
    # yMatvec = parse.(BigFloat, split(yMatstring[2:end], ","))

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
        line = nextline(io)
        line = nextline(io)
        line = nextline(io)
        structurevec = parse.(Int, split(line))
    end
    xMatbm = SDOI.BlockMatrix{BigFloat}(map(n -> zeros(BigFloat, abs(n), abs(n)), optimizer.blkdims))
    yMatbm = SDOI.BlockMatrix{BigFloat}(map(n -> zeros(BigFloat, abs(n), abs(n)), optimizer.blkdims))
    for i in 1:length(structurevec)
        dim = structurevec[i]
        if dim < 0
            dim = abs(dim)
            xblock = xMatvec[end-dim+1:end]
            deleteat!(xMatvec, length(xblock)-dim+1:length(xblock))
            xblock = Diagonal(xblock)
            xMatbm.blocks[i] = xblock

            yblock = yMatvec[end - dim + 1:end]
            deleteat!(yMatvec, length(yblock)-dim+1:length(yblock))
            yblock = Diagonal(yblock)
            yMatbm.blocks[i] = yblock


        elseif dim > 0
            xblock = xMatvec[end - dim^2 + 1:end]
            deleteat!(xMatvec, length(xblock) - dim^2 + 1:length(xblock))
            xblock = reshape(xblock, dim, dim)
            xMatbm.blocks[i] = xblock

            yblock = yMatvec[end - dim^2 + 1:end]
            deleteat!(yMatvec, length(yblock)-dim^2 +1:length(yblock))
            yblock = reshape(yblock, dim, dim)
            yMatbm.blocks[i] = yblock

        end
    end
    optimizer.X = xMatbm
    optimizer.Z = yMatbm

end
