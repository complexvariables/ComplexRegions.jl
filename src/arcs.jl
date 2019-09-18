# Type  
"""
	(type) Arc{T<:AnyComplex} in the complex plane 

Each `Arc` type is parameterized according to the common type of its complex input arguments. 
"""
struct Arc{T<:AnyComplex} <: AbstractCurve 
	circle::Circle{T} 
	start::Float64  # specified as positive fraction of 1 ccw rotation from positive real
	delta::Float64 
end

# Untyped constructors
"""
	Arc(C,start,delta)

Consruct an arc that is the part of the Circle `C` starting at parameter value `start` and ending at `start+delta`. The values are expressed in terms of fractions of a complete circle. The `start` value should be in [0,1), and `delta` should be in [-1,1].

	Arc(a,b,c)

Construct the arc starting at point `a`, passing through `b`, and ending at `c`. If the three points are collinear, a `Segment` is returned instead.
"""
function Arc(C::Circle{T},start::Real,delta::Real) where T<:AnyComplex
	if delta < 0
		Arc{T}(reverse(C),float(-start),float(-delta))
	else
		Arc{T}(C,float(start),float(delta))
	end
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
		if !C.ccw 
			ti = 1-ti
			delta = 1-delta 
		end	
		Arc(C,ti,delta)
	end
end

# Complex type converters
for ctype in [:Spherical,:Polar,:Complex]
	@eval begin
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
	t = mod(tc-A.start+DEFAULT[:tol],1)-DEFAULT[:tol]
	t/A.delta
end

unittangent(A::Arc,t::Real) = unittangent(A.circle,A.start+t*A.delta)
function tangent(A::Arc{T},t::Real) where T <: AnyComplex
	T( tangent(A.circle,A.start+t*A.delta)*A.delta )
end

+(A::Arc,z::Number) = Arc(A.circle+z,A.start,A.delta)
function -(A::Arc)
	ti = mod(A.start+0.5,1)
	Arc(-A.circle,ti,A.delta)
end

function *(A::Arc,z::Number)
	phi = angle(z)/(2*pi)
	ti = A.circle.ccw ? mod(A.start+phi,1) : mod(A.start-phi,1)
	Arc(z*A.circle,ti,A.delta)
end

"""
	inv(A) 

Invert the arc `A` through the origin. In general the inverse is an `Arc`, though the result is a `Segment` if the arc's circle passes through the origin.
"""
function inv(A::Arc) 
	w = 1 ./ point(A,[0,0.5,1])
	Arc(w...)
end

# Other methods
isfinite(::Arc) = true
conj(A::Arc) = Arc(conj(A(0)),conj(A(0.5)),conj(A(1)))
reverse(A::Arc) = Arc(A(1),A(0.5),A(0))

"""
	isapprox(A1::Arc,A2::Arc; tol=<default>) 
	A1 ≈ A2 

Determine if `A1` and `A2` represent the same arc, irrespective of the type or values of its parameters. Identity is determined by agreement within `tol`, which is interpreted as the weaker of absolute and relative differences.
"""
function isapprox(A1::Arc,A2::Arc;tol=DEFAULT[:tol])
	if isapprox(A1.circle,A2.circle,tol=tol) 
		z1 = point(A1,[0,1])
		z2 = point(A2,[0,1])
		appx(u,v) = isapprox(u,v,atol=tol,rtol=tol)
		return appx(z1,z2) || appx(z1,reverse(z2))
	else
		return false
	end
end

""" 
	dist(z,A::Arc) 

Compute the distance from number `z` to the arc `A`. 
"""
function dist(z::Number,A::Arc) 
	C = A.circle
	# result depends on whether the radius to z intersects the arc
	ti,del = A.start,A.delta
	if mod(arg(C,z)-ti,1) ≤ del 
		return dist(z,C)
	else
		return min( abs(z-point(A,0)), abs(z-point(A,1)) )
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
