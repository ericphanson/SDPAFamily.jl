# Inspired from `MOI.Utilities.ModelFilter`

struct ObjectiveFunctionFilter{T} <: MOI.ModelLike
    inner::T
    function ObjectiveFunctionFilter(model::MOI.ModelLike)
        return new{typeof(model)}(model)
    end
end

function MOI.get(model::ObjectiveFunctionFilter, attr::MOI.ObjectiveFunction)
    func = MOI.get(model.inner, attr)
    constant = MOI.constant(func)
    if !iszero(constant)
        func = copy(func)
        func.constant = zero(constant)
    end
    return MOI.Utilities.canonical(func)
end

# These just forward the attributes into the inner model.

function MOI.get(model::ObjectiveFunctionFilter, attr::MOI.AbstractModelAttribute)
    return MOI.get(model.inner, attr)
end

function MOI.get(
    model::ObjectiveFunctionFilter,
    attr::MOI.AbstractVariableAttribute,
    x::MOI.VariableIndex,
)
    return MOI.get(model.inner, attr, x)
end

function MOI.get(
    model::ObjectiveFunctionFilter,
    attr::MOI.AbstractConstraintAttribute,
    ci::MOI.ConstraintIndex,
)
    return MOI.get(model.inner, attr, ci)
end
