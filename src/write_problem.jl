using MathOptInterface
using SemidefiniteOptInterface
using Convex
const MOI = MathOptInterface
const MOIU = MOI.Utilities
const MOIB = MOI.Bridges
const SDOI = SemidefiniteOptInterface

export write_problem

MOIU.@model(MOIModel,
            (MOI.ZeroOne, MOI.Integer),
            (),
            (MOI.Zeros, MOI.Nonnegatives, MOI.Nonpositives, MOI.SecondOrderCone,
             MOI.RotatedSecondOrderCone,
             MOI.PositiveSemidefiniteConeTriangle,),
            (),
            (MOI.SingleVariable,),
            (MOI.ScalarAffineFunction,),
            (MOI.VectorOfVariables,),
            (MOI.VectorAffineFunction,))

# function load_problem_mock!(problem::Problem{T}, optimizer::MOI.ModelLike) where {T}
#
#     model = MOIU.CachingOptimizer(MOIModel{T}(), MOIU.MANUAL)
#     var_to_ranges = Convex.load_MOI_model!(model, problem)
#
#     universal_fallback = MOIU.UniversalFallback(MOIModel{T}())
#     optimizer = MOIU.CachingOptimizer(universal_fallback, optimizer)
#     optimizer = MOI.Bridges.full_bridge_optimizer(optimizer, T)
#
#     MOIU.reset_optimizer(model, optimizer);
#     MOIU.attach_optimizer(model);
#     MOI.optimize!(model)
#
#     return optimizer # not necessary?
# end
"""
    function write_problem(problem::Problem{T}, filename::String)
Write the _dual_ problem parameters into the specified file. The file path is returned.

"""
function write_problem(mock::SDOI.MockSDOptimizer{T}, problem::Problem{T}, filepath::String) where {T}
    # define mock optimizer
    mock_optimizer = SDOI.SDOIOptimizer(mock, T)

    # lines borrowed from Convex.solve!() to load problem into mock_optimizer
    model = MOIU.CachingOptimizer(MOIModel{T}(), MOIU.MANUAL)
    var_to_ranges = Convex.load_MOI_model!(model, problem)

    universal_fallback = MOIU.UniversalFallback(MOIModel{T}())
    optimizer = MOIU.CachingOptimizer(universal_fallback, mock_optimizer)
    optimizer = MOI.Bridges.full_bridge_optimizer(optimizer, T)

    MOIU.reset_optimizer(model, optimizer);
    MOIU.attach_optimizer(model);
    MOI.optimize!(model)

    MOI.write(mock_optimizer, filepath)
    # relativepath = joinpath("./", filename)
    # filepath = abspath(relativepath)
    # @info "Problem written to file $filepath"
    return filepath
end

# function write_problem_trycatch(problem::Problem{T}, filename::String) where {T}
#
#     # another hacky way to write the problem...
#
#     mock = SDOI.MockSDOptimizer{T}()
#     mock_optimizer = SDOI.SDOIOptimizer(mock, T)
#     try
#         Convex.solve!(problem, mock_optimizer)
#     catch
#     end
#     MOI.write(mock_optimizer, filename)
# end
