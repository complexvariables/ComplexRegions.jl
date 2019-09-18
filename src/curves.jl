#
# abstract interfaces
#

abstract type AbstractCurve end

# Required methods
"""
	point(C::AbstractCurve,t::Real)
Find the point on curve `C` at parameter value `t`, which should lie in the interval [0,1]. 
"""
point(c::AbstractCurve,t::Real) = @error "No point() method defined for type $(typeof(c))"

"""
	tangent(C::AbstractCurve,t::Real)
Find the complex number representing the tangent to curve `C` at parameter value `t` in [0,1]. 
"""
tangent(C::AbstractCurve,t::Real) = @error "No tangent() method defined for type $(typeof(C))"

reverse(C::AbstractCurve) = @error "No reverse() method defined for type $(typeof(C))"

""" 
	isfinite(C::AbstractCurve) 
Return `true` if the curve is bounded in the complex plane (i.e., does not pass through infinity).
"""
isfinite(C::AbstractCurve) = @error "No isfinite() method defined for type $(typeof(C))"

conj(C::AbstractCurve) = @error "No conj() method defined for type $(typeof(C))"
+(C::AbstractCurve,z::Number) = @error "No addition method defined for type $(typeof(C))"
-(C::AbstractCurve) = @error "No negation method defined for type $(typeof(C))"
*(C::AbstractCurve,z::Number) = @error "No multiplication method defined for type $(typeof(C))"
inv(C::AbstractCurve) = @error "No inversion method defined for type $(typeof(C))"

# Default implementations
"""
	point(C::AbstractCurve,t::AbstractArray)

Vectorize the `point` function for curve `C`.
"""
point(c::AbstractCurve,t::AbstractArray{T}) where T<:Real = [point(c,t) for t in t]

"""
	unittangent(C::AbstractCurve,t::Real)
Find the complex number representing the unit tangent to curve `C` at parameter value `t` in [0,1]. For Lines, Segments, and Rays, the `t` argument is optional.
"""
unittangent(C::AbstractCurve,t::Real) = sign(tangent(C,t))

"""
	normal(C::AbstractCurve,t::Real)
Find the unit complex number in the direction of the leftward-pointing normal to curve `C` at parameter value `t` in [0,1]. 
"""
normal(c::AbstractCurve,t::Real) = 1im*unittangent(c,t)

+(z::Number,C::AbstractCurve) = +(C,z)
-(C::AbstractCurve,z::Number) = C + (-z)
-(z::Number,C::AbstractCurve) = z + (-C)
*(z::Number,C::AbstractCurve) = *(C,z)
/(C::AbstractCurve,z::Number) = C*(1/z)
/(z::Number,C::AbstractCurve) = z*inv(C)

function arclength(C::AbstractCurve,part=[0,1])
	f = t -> abs(tangent(C,t))
	intadapt(f,part...,DEFAULT[:tol])
end

"""
	plotdata(C::AbstractCurve,n=501)

Compute `n` points along the curve `C` suitable to make a plot of it.
"""
plotdata(C::AbstractCurve) = adaptpoints(t->point(C,t),t->unittangent(C,t),0,1)

show(io::IO,C::AbstractCurve) = print(io,"Complex-valued $(typeof(C))")
show(io::IO,::MIME"text/plain",C::AbstractCurve) = print(io,"Complex-valued $(typeof(C))")

# AbstractClosedCurve
abstract type AbstractClosedCurve <: AbstractCurve end

# Default implementations
function winding(C::AbstractClosedCurve,z::Number)
	# Integrate around the curve
	f = t -> imag(tangent(C,t)/(point(C,t)-z))
	w = intadapt(f,0,1,1e-4)
	return round(Int,w/(2π))
end

isinside(z::Number,C::AbstractClosedCurve) = winding(C,z) != 0
isoutside(z::Number,C::AbstractClosedCurve) = winding(C,z) == 0

#
# generic Curve 
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

	Curve(f,df[,a,b])
Construct a curve with point location and tangent given by the complex-valued functions `f` and `df`, respectively, optionally with given limits on the parameter.
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
isfinite(C::Curve) = true

+(C::Curve,z::Number) = Curve(t->C.point(t)+z,C.tangent)
-(C::Curve) = Curve(t->-C.point(t),t->-C.tangent(t))
*(C::Curve,z::Number) = Curve(t->C.point(t)*z,t->C.tangent(t)*z)
inv(C::Curve) = Curve(t->1/C.point(t),t->-C.tangent(t)/C.point(t)^2)

#
# generic ClosedCurve
# 

"""
(type) Smooth closed curve defined by an explicit function of a real paramerter in [0,1]. 
"""
struct ClosedCurve <: AbstractClosedCurve 
	curve::Curve
	function ClosedCurve(c::Curve;tol=DEFAULT[:tol])
		@assert isapprox(point(c,0),point(c,1);rtol=tol,atol=tol) "Curve does not close"
		new(c)
	end
end

"""
	ClosedCurve(f; tol=<default>)
	ClosedCurve(f,a,b; tol=<default)
Construct a `ClosedCurve` object from the complex-valued function `point` accepting an argument in the interval [0,1]. The constructor checks whether `f(0)≈f(1)` to tolerance `tol`. If `a` and `b` are given, they are the limits of the parameter in the call to the supplied `f`. However, the resulting object will be defined on [0,1], which is internally scaled to [a,b].

	ClosedCurve(f,df[,a,b]; tol=<default>)
Construct a closed curve with point location and tangent given by the complex-valued functions `f` and `df`, respectively, optionally with given limits on the parameter.
"""
ClosedCurve(f,df=t->fdtangent(f,t);kw...) = ClosedCurve(Curve(f,df;kw...))
ClosedCurve(f,a::Real,b::Real;kw...) = ClosedCurve(Curve(f,a,b;kw...))
ClosedCurve(f,df,a::Real,b::Real;kw...) = ClosedCurve(Curve(f,df,a,b;kw...))

point(C::ClosedCurve,t::Real) = point(C.curve,t)
(C::ClosedCurve)(t::Real) = point(C.curve,t)
tangent(C::ClosedCurve,t::Real) = tangent(C.curve,t)
conj(C::ClosedCurve) = ClosedCurve(conj(C.curve))
reverse(C::ClosedCurve) = ClosedCurve(reverse(C.curve))
isfinite(C::ClosedCurve) = isfinite(C.curve) 
+(C::ClosedCurve,z::Number) = ClosedCurve(C.curve+z)
-(C::ClosedCurve) = ClosedCurve(-C.curve)
*(C::ClosedCurve,z::Number) = ClosedCurve(C.curve*z)
inv(C::ClosedCurve) = ClosedCurve(inv(C.curve))

include("lines.jl")
include("rays.jl")
include("segments.jl")
include("circles.jl")
include("arcs.jl")
include("intersections.jl")
