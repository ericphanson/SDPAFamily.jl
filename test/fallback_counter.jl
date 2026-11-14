using Logging

# A logger that counts the fallbacks of the presolve (debug messages with the group
# `:presolve_fallback`) and passes all other messages to `parent`.
mutable struct FallbackCounter <: Logging.AbstractLogger
    parent::Logging.AbstractLogger
    reasons::Vector{String}
end
FallbackCounter(parent = Logging.current_logger()) = FallbackCounter(parent, String[])

Logging.min_enabled_level(::FallbackCounter) = Logging.Debug
function Logging.shouldlog(l::FallbackCounter, level, _module, group, id)
    group === :presolve_fallback && return true
    return level >= Logging.min_enabled_level(l.parent) &&
           Logging.shouldlog(l.parent, level, _module, group, id)
end
function Logging.handle_message(l::FallbackCounter, level, message, _module, group, id,
                                file, line; kwargs...)
    if group === :presolve_fallback
        push!(l.reasons, string(get(kwargs, :reason, "")))
        return nothing
    end
    return Logging.handle_message(l.parent, level, message, _module, group, id, file, line; kwargs...)
end
Logging.catch_exceptions(l::FallbackCounter) = Logging.catch_exceptions(l.parent)
