abstract type AbstractCurve end

point(c::AbstractCurve,t::Real) = @error "No point() method defined for type $(typeof(c))"
point(c::AbstractCurve,t::AbstractArray{T}) where T<:Real = [point(c,t) for t in t]
start(c::AbstractCurve) = point(c,0.0)
stop(c::AbstractCurve) = point(c,1.0)
arclength(c::AbstractCurve) = @error "No arclength() method defined for type $(typeof(c))"

struct Curve <: AbstractCurve 
	point 
	arclength 
end
Curve(f) = Curve(f,missing)
Curve(f,a::Real,b::Real,arclen=missing) = Curve(t -> f(scaleto(a,b,t)),arclen)

point(C::Curve,t::Real) = C.point(t)
(C::Curve)(t::Real) = point(C,t)
arclength(C::Curve) = C.arclength 

abstract type AbstractClosedCurve <: AbstractCurve end

struct ClosedCurve <: AbstractClosedCurve 
	point 
	arclength 
	function ClosedCurve(f,arclen=missing;tol=1e-12)
		@assert isapprox(f(0),f(1);rtol=tol,atol=tol) "Curve does not close"
		new(f,arclen)
	end
end
ClosedCurve(f,a::Real,b::Real,arclen=missing;kw...) = ClosedCurve(t -> f(scaleto(a,b,t)),arclen;kw...)

point(C::ClosedCurve,t::Real) = C.point(t)
(C::ClosedCurve)(t::Real) = point(C,t)
arclength(C::ClosedCurve) = C.arclength 

include("lines_segments.jl")
include("circles_arcs.jl")
