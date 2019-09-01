# Type and constructors
"""
	(type) Circle{T<:AnyComplex} in the complex plane 

Each `Circle` type is parameterized according to the common type of its complex input arguments. 
"""
struct Circle{T<:AnyComplex} <: AbstractClosedCurve 
	center::T
	radius::Float64
	ccw::Bool
end
"""
	Circle(zc,r,ccw=true)

Construct the circle with given center `zc`, radius `r`, and orientation (defaults to counterclockwise).

	Circle(a,b,c)

Construct the circle passing through the given three numbers. Orientation is determined so that the values are visited in the given order. If the three points are collinear (including when one of the given values is infinite), a `Line` is returned instead.
"""
Circle(z::AnyComplex,r::Real,ccw::Bool=true) = Circle{typeof(z)}(z,r,ccw)
Circle(z::Number,r::Real,ccw::Bool=true) = Circle(complex(float(z)),r,ccw)

# Construction by three points
function Circle(a::Number,b::Number,c::Number) 
	a,b,c = promote( complex.(float.([a,b,c]))... )
	Circle(a,b,c)
end
function Circle(a::T,b::T,c::T) where {T<:AnyComplex}
	isinf(a) && return Line(b,c)
	isinf(b) && return Line(c,a)
	isinf(c) && return Line(a,b)
	# Use intersection of chord bisectors to find the center of the circle. 
	w = (a-c)/2
	d1,d2 = a-b,c-b
	M = SMatrix{2,2}(real(d1),imag(d1),real(d2),imag(d2))
	if cond(M) > 0.1/eps(typeof(float(real(a))))
		# Collinear points
		return Line(a,b)
	else
		p =  M \ SVector(imag(w),-real(w))
		cen = (a+b)/2 - 1im*p[1]*d1
		ccw = isccw(a-cen,b-cen,c-cen)
		return Circle{T}(cen,abs(a-cen),ccw)
	end
end

# Complex type converters
for ctype in [:Spherical,:Polar,:Complex]
	@eval begin
		function $ctype(C::Circle{T}) where T<:AnyComplex 
			Circle($ctype(C.center),C.radius,C.ccw)
		end	
	end
end

# Required methods
function point(C::Circle,t::Real)
	z = C.ccw ? exp(2im*pi*t) : exp(-2im*pi*t)
	C.center + typeof(C.center)(z*C.radius)
end
(C::Circle)(t::Real) = point(C,t)

ispositive(C::Circle) = C.ccw

arclength(C::Circle) = 2π*C.radius

""" 
	arg(C::Circle,z) 

Find the parameter argument `t` such that `C(t)==z` is true. 

This gives undefined results if `z` is not actually on the circle. 
"""
function arg(C::Circle,z::Number)
	α = angle(z-C.center)/(2π)
	C.ccw ? mod(α,1) : mod(-α,1)
end

function unittangent(C::Circle{T},t::Real) where T <: AnyComplex
	C.ccw ? T( 1im*exp(2im*pi*t) ) : T( -1im*exp(-2im*pi*t) )
end
function tangent(C::Circle{T},t::Real) where T <: AnyComplex
	T( 2π*C.radius*unittangent(C,t) )
end

+(C::Circle,z::Number) = Circle(C.center+z,C.radius,C.ccw)
-(C::Circle) = Circle(-C.center,C.radius,C.ccw)
*(C::Circle,z::Number) = Circle(C.center*z,C.radius*abs(z),C.ccw)

"""
	inv(C) 
Invert the circle `C` through the origin. In general the inverse is a `Circle`, though the result is a `Line` if `C` passes through the origin.
"""
function inv(C::Circle) 
	w = 1 ./ point(C,[0,0.25,0.5])
	Circle(w...)
end

# Other methods
isfinite(::Circle) = true 
conj(C::Circle) = Circle(conj(C.center),C.radius,!C.ccw)
reverse(C::Circle) = Circle(C.center,C.radius,!C.ccw)

"""
	isapprox(C1::Circle,C2::Circle; tol=<default>) 
	C1 ≈ C2 
Determine if `C1` and `C2` represent the same circle, irrespective of the type or values of its parameters. Identity is determined by agreement within `tol`, which is interpreted as the weaker of absolute and relative differences.
"""
function isapprox(C1::Circle,C2::Circle;tol=DEFAULT[:tol])
	return isapprox(C1.center,C2.center,rtol=tol,atol=tol) &&
	isapprox(C1.radius,C2.radius,rtol=tol,atol=tol)
end

function winding(C::Circle,z::Number)
	w =  abs(z-C.center) < C.radius ? 1 : 0
	C.ccw ? w : -w
end

""" 
	dist(z,C::Circle) 
Compute the distance from number `z` to the circle `C`. 
"""
dist(z::Number,C::Circle) = abs( abs(z-C.center) - C.radius )

""" 
	closest(z,C::Circle) 
Find the point on circle `C` that lies closest to `z`.
"""
closest(z::Number,C::Circle) =	C.center + C.radius*sign(z - C.center)

""" 
	reflect(z,C::Circle) 
Reflect the value `z` across the circle `C`. (For reflection of a circle through a point, use translation and negation.)
"""
function reflect(z::Number,C::Circle)
	ζ = z-C.center
	ζ==0 ? convert(typeof(float(z)),Inf) : C.center + C.radius^2/conj(ζ)
end

function show(io::IO,C::Circle)
	orient = C.ccw ? "ccw" : "cw"
	print(IOContext(io,:compact=>true),"Circle(",C.center,",",C.radius,",",orient,")")
end
function show(io::IO,::MIME"text/plain",C::Circle{T}) where {T}
	orient = C.ccw ? "positively" : "negatively"
	print(io,"Circle{$T} in the complex plane:\n   centered at (",C.center,") with radius $(C.radius), $(orient) oriented")
end
