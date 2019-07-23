# Type  
"""
	(type) Arc{T<:AnyComplex} in the complex plane 

Each `Arc` type is parameterized according to the common type of its input arguments. 
"""
struct Arc{T<:AnyComplex} <: AbstractCurve 
	circle::Circle{T} 
	start::Float64  # specified as positive fraction of 1 ccw rotation from positive real
	delta::Float64 
end

# Untyped constructors
"""
	Arc(c,start,delta)

Consruct an arc that is the part of the Circle `c` starting at parameter value `start` and ending at `start+delta`. The values are expressed in terms of fractions of a complete circle. The `start` value should be in [0,1), and `delta` should be in [-1,1].

	Arc(a,b,c)

Construct the arc starting at point `a`, passing through `b`, and ending at `c`. If the three points are collinear, a `Segment` is returned instead.

	Arc(a,b,center=zc)

Construct the arc starting at point `a`, ending at `b`, and whose circular center is at `zc`. No check is performed to see if `b` actually lies on the circle. 
"""
function Arc(C::Circle{T},start::Real,delta::Real) where T<:AnyComplex
	Arc{T}(C,Float64(start),Float64(delta))
end

# Construct from 3 points
function Arc(a::Number,m::Number,b::Number) 
	a,m,b = promote(complex(float(a)),m,b)
	C = Circle(a,m,b)
	if isa(C,Line)  # collinear
		Segment(a,b)
	else
		α,β = a-C.center,b-C.center
		ti = mod(angle(α)/(2π),1)
		delta = mod(angle(β/α)/(2π),1) # force into (0,1)
		# which of the two circle pieces do we use? 
		if mod(angle((m-C.center)/α)/(2π),1) > delta 
			delta = delta-1
		end	
		Arc(C,ti,delta)
	end
end

# Construct from 2 points and circle center 
function Arc(a::Number,b::Number;center=0) 
	a,b,zc = promote(complex(float(a)),b,center)
	C = Circle(zc,abs(a-zc))
	if isa(C,Line)  # collinear
		Segment(a,b)
	else
		α = a-C.center
		ti = mod(angle(α)/(2π),1)
		delta = angle((b-C.center)/α)/(2π)
		Arc(C,ti,delta)
	end
end

# Complex type converters
for ctype in [:Spherical,:Polar,:Complex]
	docstr = """
	$ctype(::Arc)
Convert to `Arc{$ctype}`. This is useful for plotting curves in a desired way.
"""
	@eval begin
		@doc $docstr -> 
		function $ctype(A::Arc{T}) where T<:AnyComplex 
			Arc($ctype(A.circle),A.start,A.delta)
		end	
	end
end

# Required methods
function point(A::Arc,t::Real) 
	s = scaleto(A.start,A.start+A.delta,t)
	point(A.circle,s)
end
arclength(A::Arc) = arclength(A.circle)*A.delta
(C::Arc)(t::Real) = point(C,t)
""" 
	arg(A::Arc,z) 

Find the parameter argument `t` such that `A(t)==z` is true. 

This gives undefined results if `z` is not actually on the arc. 
"""
function arg(A::Arc,z::Number)
	tc = arg(A.circle,z)
	t = mod(tc-A.start,1)
	A.delta < 0 ? -mod(1-t,1)/A.delta : t/A.delta
end
tangent(A::Arc,t::Real) = tangent(A.circle,A.start + t*A.delta)

# Other methods
isfinite(::Arc) = true
conj(A::Arc) = Arc(conj(A(0)),conj(A(0.5)),conj(A(1)))
reverse(A::Arc) = Arc(A(1),A(0.5),A(0))
"""
	A + z
	z + A 

Translate the arc `A` by a number `z`. 
"""
+(A::Arc,z::Number) = Arc(A.circle+z,A.start,A.delta)
+(z::Number,A::Arc) = Arc(z+A.circle,A.start,A.delta)
"""
	A - z

Translate the arc `A` by a number `-z`.

	-A 
	z - A 

Negate a arc `A` (reflect through the origin), and optionally translate by a number `z`.
"""
-(A::Arc,z::Number) = Arc(A.circle-z,A.start,A.delta)
function -(A::Arc)
	ti = mod(A.start+0.5,1)
	Arc(-A.circle,ti,A.delta)
end
-(z::Number,A::Arc) = z + (-A)
"""
	z*A 
	A*z 

Multiply the arc `A` by real or complex number `z`; i.e., scale and rotate it about the origin.
"""
function *(A::Arc,z::Number)
	phi = angle(z)/(2*pi)
	ti = mod(A.start+phi,1)
	Arc(z*A.circle,ti,A.delta)
end
*(z::Number,A::Arc) = A*z
"""
	A/z 

Multiply the arc `A` by the number `1/z`; i.e., scale and rotate it about the origin.

	z/A 
	inv(A) 

Invert the arc `A` through the origin (and optionally multiply by the number `1/z`). In general the inverse is an `Arc`, though the result is a `Segment` if the arc's circle passes through the origin.
"""
/(A::Arc,z::Number) = A*(1/z)
function /(z::Number,A::Arc) 
	w = z./point(A,[0,0.5,1])
	Arc(w...)
end
inv(A::Arc) = 1/A

"""
	isapprox(A1::Arc,A2::Arc; tol=<default>) 
	A1 ≈ A2 

Determine if `A1` and `A2` represent the same arc, irrespective of the type or values of its parameters. Identity is determined by agreement within `tol`, which is interpreted as the weaker of absolute and relative differences.
"""
function isapprox(A1::Arc,A2::Arc;tol=DEFAULT[:tol])
	return isapprox(A1.Circle,A2.Circle,tol) &&
		isapprox(A1.start,A2.start,rtol=tol,atol=tol) &&
		isapprox(A1.delta,A2.delta,rtol=tol,atol=tol) 
end

""" 
	dist(z,A::Arc) 

Compute the distance from number `z` to the arc `A`. 
"""
function dist(z::Number,A::Arc) 
	if A.delta > 0
		ti,del = A.start,A.delta
	else
		ti = mod(A.start+A.delta,1)
		del = -A.delta 
	end
	ζ = z - A.circle.center
	α = mod(angle(ζ)/(2π)-ti,1)
	if 0 ≤ α ≤ del 
		return abs(abs(ζ)-A.circle.radius)
	else
		return min(abs(z-point(A,0)), abs(z-point(A,1)) )
	end
end
""" 
	closest(z,A::Arc) 

Find the point on arc `A` that lies closest to `z`.
"""
function closest(z::Number,A::Arc)
	ζ = z - A.circle.center
	d = A.delta/2
	# rotate arc to position symmetric about positive Re axis
	ϕ = angle( ζ/exp(2im*pi*(A.start+d)) ) / (2π)
	if ϕ > d 
		point(A,1)
	elseif ϕ < -d 
		point(A,0) 
	else
		A.circle.center + A.circle.radius*sign(ζ)
	end
end

function show(io::IO,A::Arc{T}) where {T}
	print(IOContext(io,:compact=>true),"Arc(",A(0.0),"...",A(1.0),")")
end
function show(io::IO,::MIME"text/plain",A::Arc{T}) where {T}
	print(io,"Arc{$T} in the complex plane:\n   fraction ",A.delta," of (",A.circle,") starting at ",A.start)
end
