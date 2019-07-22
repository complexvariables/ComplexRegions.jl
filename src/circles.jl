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
point(C::Circle,t::Real) = C.center + C.radius*exp(2im*pi*t)
arclength(C::Circle) = 2π*C.radius
(C::Circle)(t::Real) = point(C,t)
arg(C::Circle,z::Number) = mod(angle(z-C.center)/(2π),1)
tangent(C::Circle,t::Real) = C.ccw ? 1im*exp(2im*pi*t) : -1im*exp(2im*pi*t)

# Other methods
isfinite(::Circle) = true 
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
inv(C::Circle) = 1/C

isleft(z::Number,C::Circle) = !xor(C.ccw,abs(z-C.center) < C.radius) 

function isapprox(C1::Circle,C2::Circle;tol=DEFAULT[:tol])
	return isapprox(C1.center,C2.center,rtol=tol,atol=tol) &&
		isapprox(C1.radius,C2.radius,rtol=tol,atol=tol)
end

dist(z::Number,C::Circle) = abs( abs(z-C.center) - C.radius )
closest(z::Number,C::Circle) =	C.center + C.radius*sign(z - C.center)
function reflect(z::Number,C::Circle)
	ζ = z-C.center
	ζ==0 ? convert(typeof(z),Inf) : C.center + ζ/abs2(ζ)
end

function show(io::IO,C::Circle)
	print(IOContext(io,:compact=>true),"Circle(",C.center,",",C.radius,")")
end
function show(io::IO,::MIME"text/plain",C::Circle{T}) where {T}
	print(io,"Circle{$T} in the complex plane:\n   centered at (",C.center,") with radius ",C.radius)
end
