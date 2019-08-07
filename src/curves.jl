#
# abstract interface
#

abstract type AbstractCurve end
abstract type AbstractClosedCurve <: AbstractCurve end

"""
	point(C::AbstractCurve,t::Real)

Find the point on curve `C` at parameter value `t`, which should lie in the interval [0,1]. 
"""
point(c::AbstractCurve,t::Real) = @error "No point() method defined for type $(typeof(c))"

"""
	point(C::AbstractCurve,t::AbstractArray)

Vectorize the `point` function for curve `C`.
"""
point(c::AbstractCurve,t::AbstractArray{T}) where T<:Real = [point(c,t) for t in t]
# """
# 	arclength(C::AbstractCurve)

# Fetch or compute the arc length of curve `C`.
# """
# arclength(c::AbstractCurve) = @error "No arclength() method defined for type $(typeof(c))"
"""
	tangent(C::AbstractCurve,t::Real)

Find the complex number representing the tangent to curve `C` at parameter value `t` in [0,1]. 
"""
tangent(C::AbstractCurve,t::Real) = @error "No tangent() method defined for type $(typeof(C))"

"""
	unittangent(C::AbstractCurve,t::Real)

Find the complex number representing the unit tangent to curve `C` at parameter value `t` in [0,1]. For Lines, Segments, and Rays, the `t` argument is optional.
"""
unittangent(C::AbstractCurve,t::Real) = sign(tangent(C,t))

"""
	normal(C::AbstractCurve,t::Real)

Find the unit complex number in the direction of the leftward-pointing normal to curve `C` at parameter value `t` in [0,1]. 
"""
normal(c::AbstractCurve,t::Real) = 1im*tangent(c,t)

"""
	plotdata(C::AbstractCurve,n=501)

Compute `n` points along the curve `C` suitable to make a plot of it.
"""
plotdata(C::AbstractCurve) = adaptpoints(t->point(C,t),t->unittangent(C,t),0,1)

# """
# 	conj(C::AbstractCurve) 

# Construct the curve conjugate to `C`. (For a closed curve, the result reverses the orientation with respect to its interior.)
# """
# conj(C::AbstractCurve) = @error "No conj() method defined for type $(typeof(C))"
"""
	reverse(C::AbstractCurve)

Construct a curve identical to `C` except for the direction of traversal. (Essentially, replace "t" by "1-t" in the parameterization.)
"""
reverse(C::AbstractCurve) = @error "No reverse() method defined for type $(typeof(C))"

""" 
	isfinite(C::AbstractCurve) 

Return `true` if the curve is bounded in the complex plane (i.e., does not pass through infinity).
"""
isfinite(C::AbstractCurve) = @error "No isfinite() method defined for type $(typeof(C))"

#
# generic curve type 
#

"""
(type) Smooth curve defined by an explicit function of a real paramerter in [0,1]. 
"""
struct Curve <: AbstractCurve 
	point 
	tangent
end

"""
	Curve(f)
	Curve(f,a,b)

Construct a `Curve` object from the complex-valued function `f` accepting an argument in the interval [0,1]. If `a` and `b` are given, they are the limits of the parameter in the call to the supplied `f`. However, the resulting object will be defined on [0,1], which is internally scaled to [a,b].
"""
Curve(f) = Curve(f,t->fdtangent(f,t))
Curve(f,df,a::Real,b::Real) = Curve(t->f(scaleto(a,b,t)),t->df(scaleto(a,b,t)))
Curve(f,a::Real,b::Real) = Curve(t->f(scaleto(a,b,t)))

# Required methods
point(C::Curve,t::Real) = C.point(t)
(C::Curve)(t::Real) = point(C,t)
tangent(C::Curve,t::Real) = C.tangent(t)
conj(C::Curve) = Curve(t->conj(C.point(t)))
reverse(C::Curve) = Curve(t->C.point(1-t))

#
# generic closed curve type
# 

"""
(type) Smooth closed curve defined by an explicit function of a real paramerter in [0,1]. 
"""
struct ClosedCurve <: AbstractClosedCurve 
	point 
	tangent
	function ClosedCurve(f,df=t->fdtangent(f,t);tol=DEFAULT[:tol])
		@assert isapprox(f(0),f(1);rtol=tol,atol=tol) "Curve does not close"
		new(f,df)
	end
end

"""
	ClosedCurve(f,arclen=missing;tol=DEFAULT[:tol])
	ClosedCurve(f,a,b,arclen=missing;tol=DEFAULT[:tol])

Construct a `ClosedCurve` object from the complex-valued function `point` accepting an argument in the interval [0,1]. The constructor checks whether `f(0)≈f(1)` to tolerance `tol`. 

If `a` and `b` are given, they are the limits of the parameter in the call to the supplied `f`. However, the resulting object will be defined on [0,1], which is internally scaled to [a,b].
"""
ClosedCurve(f,a::Real,b::Real;kw...) = ClosedCurve(t -> f(scaleto(a,b,t));kw...)
function ClosedCurve(f,df,a::Real,b::Real;kw...)
	ClosedCurve(t->f(scaleto(a,b,t)), t->df(scaleto(a,b,t)); kw...)
end

point(C::ClosedCurve,t::Real) = C.point(t)
(C::ClosedCurve)(t::Real) = point(C,t)
tangent(C::ClosedCurve,t::Real) = C.tangent(t)
conj(C::ClosedCurve) = ClosedCurve(t->conj(C.point(t)))
reverse(C::ClosedCurve) = ClosedCurve(t->C.point(1-t))

include("lines.jl")
include("rays.jl")
include("segments.jl")
include("circles.jl")
include("arcs.jl")
include("intersections.jl")
