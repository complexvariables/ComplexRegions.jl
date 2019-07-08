# Circle 
# Type and constructors
struct Circle{T<:AnyComplex} <: AbstractClosedCurve 
	center::T
	radius::Float64
	ccw::Bool
end
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
	M = [ real(d1) real(d2); imag(d1) imag(d2) ]
	if cond(M) > 0.1/eps(typeof(float(real(a))))
		# Collinear points
		return Line(a,b)
	else
		p =  M \ [imag(w);-real(w)] 
		cen = (a+b)/2 - 1im*p[1]*d1
		ccw = mod2pi(angle((b-cen)/(a-cen))) < mod2pi(angle((c-cen)/(a-cen)))
		return Circle{T}(cen,abs(a-cen),ccw)
	end
end

# Required methods
point(C::Circle,t::Real) = C.center + C.radius*exp(2im*pi*t)
arclength(C::Circle) = 2π*C.radius
(C::Circle)(t::Real) = point(C,t)

# Other methods
isbounded(::Circle) = true 
conj(C::Circle) = Circle(conj(C.center),C.radius,!C.ccw)
reverse(C::Circle) = Circle(C.center,C.radius,!C.ccw)
+(C::Circle,z::Number) = Circle(C.center+z,C.radius,C.ccw)
+(z::Number,C::Circle) = Circle(C.center+z,C.radius,C.ccw)
-(C::Circle) = Circle(-C.center,C.radius,!C.ccw)
-(C::Circle,z::Number) = Circle(C.center-z,C.radius,C.ccw)
-(z::Number,C::Circle) = z + (-C)
*(C::Circle,z::Number) = Circle(C.center*z,C.radius*abs(z),C.ccw)
*(z::Number,C::Circle) = Circle(C.center*z,C.radius*abs(z),C.ccw)
/(C::Circle,z::Number) = Circle(C.center/z,C.radius/abs(z),C.ccw)
function /(z::Number,C::Circle) 
	w = z./point(C,[0,0.25,0.5])
	Circle(w...)
end

isleft(z::Number,C::Circle) = !xor(C.ccw,abs(z-C.center) < C.radius) 

function isapprox(C1::Circle,C2::Circle;tol=1e-12)
	return isapprox(C1.center,C2.center,rtol=tol,atol=tol) &&
		isapprox(C1.radius,C2.radius,rtol=tol,atol=tol)
end

dist(z::Number,C::Circle) = abs( abs(z-C.center) - C.radius )
closest(z::Number,C::Circle) =	C.center + C.radius*sign(z - C.center)

function show(io::IO,C::Circle)
	print(IOContext(io,:compact=>true),"Circle(",C.center,",",C.radius,")")
end
function show(io::IO,::MIME"text/plain",C::Circle{T}) where {T}
	print(io,"Circle{$T} in the complex plane:\n   centered at (",C.center,") with radius ",C.radius)
end

#
# Arc 
# 

# Type  
struct Arc{T<:AnyComplex} <: AbstractCurve 
	circle::Circle{T} 
	start::Float64  # specified as positive fraction of 1 ccw rotation from positive real
	delta::Float64 
end

# Untyped constructor
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
		delta = angle(β/α)/(2π)
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
		α,β = a-C.center,b-C.center
		ti = mod(angle(α)/(2π),1)
		delta = angle(β/α)/(2π)
		Arc(C,ti,delta)
	end
end

# Required methods
function point(A::Arc,t::Real) 
	s = scaleto(A.start,A.start+A.delta,t)
	point(A.circle,s)
end
arclength(A::Arc) = arclength(A.circle)*A.delta
(C::Arc)(t::Real) = point(C,t)

# Other methods
isbounded(::Arc) = true
conj(A::Arc) = Arc(conj(A(0)),conj(A(0.5)),conj(A(1)))
reverse(A::Arc) = Arc(A(1),A(0.5),A(0))
+(A::Arc,z::Number) = Arc(A.circle+z,A.start,A.delta)
+(z::Number,A::Arc) = Arc(z+A.circle,A.start,A.delta)
-(A::Arc,z::Number) = Arc(A.circle-z,A.start,A.delta)
function -(A::Arc)
	ti = mod(A.start+0.5,1)
	Arc(-A.circle,ti,A.delta)
end
-(z::Number,A::Arc) = z + (-A)
function *(A::Arc,z::Number)
	phi = angle(z)/(2*pi)
	ti = mod(A.start+phi,1)
	Arc(z*A.circle,ti,A.delta)
end
*(z::Number,A::Arc) = A*z
/(A::Arc,z::Number) = A*(1/z)
function /(z::Number,A::Arc) 
	w = z./point(A,[0,0.5,1])
	Arc(w...)
end

function isapprox(A1::Arc,A2::Arc;tol=1e-12)
	return isapprox(A1.Circle,A2.Circle,tol) &&
		isapprox(A1.start,A2.start,rtol=tol,atol=tol) &&
		isapprox(A1.delta,A2.delta,rtol=tol,atol=tol) 
end
function dist(z::Number,A::Arc) 
	ζ = z - A.circle.center 
	min( abs(abs(ζ)-A.circle.radius), abs(z-point(A,0)), abs(z-point(A,1)) )
end
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
