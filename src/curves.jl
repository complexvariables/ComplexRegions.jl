abstract type AbstractCurve end

point(c::AbstractCurve,t::Real) = @error "No point() method defined for type $(typeof(c))"
point(c::AbstractCurve,t::AbstractArray{T}) where T<:Real = [point(c,t) for t in t]
start(c::AbstractCurve) = point(c,0.0)
stop(c::AbstractCurve) = point(c,1.0)
arclength(c::AbstractCurve) = @error "No arclength() method defined for type $(typeof(c))"

abstract type AbstractClosedCurve <: AbstractCurve end

include("lines_segments.jl")
include("circles_arcs.jl")
