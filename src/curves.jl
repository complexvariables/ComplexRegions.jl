abstract type AbstractCurve end

point(c::AbstractCurve,t) = @error "No point() method defined for type $(typeof(c))"
start(c::AbstractCurve) = point(c,0.0)
stop(c::AbstractCurve) = point(c,1.0)
arclength(c::AbstractCurve) = @error "No arclength() method defined for type $(typeof(c))"

abstract type AbstractClosedCurve <: AbstractCurve end

for T in [Complex,Polar,Spherical]
	struct Circle{T} <: AbstractClosedCurve 
		center::T 
		radius::Float64
	end
end
Circle(z::Complex,r::Real) = Circle{Complex}(z,r)
Circle(z::Polar,r::Real) = Circle{Polar}(z,r)
Circle(z::Spherical,r::Real) = Circle{Spherical}(z,r)
Circle(a::Number,b::Number,c::Number) = Circle(promote(a,b,c)...)
function Circle(a::T,b::T,c::T) where {T<:AnyComplex}
	isinf(a) && return Line(b,c-b)
	isinf(b) && return Line(a,c-a)
	isinf(c) && return Line(a,b-a)
	# Use intersection of chord bisectors to find the center of the circle. 
	w = (a-c)/2
	d1,d2 = a-b,c-b
	M = [ real(d1) real(d2); imag(d1) imag(d2) ]
	if cond(M) > 0.1/eps(typeof(float(real(a))))
		# Collinear points
		return Line(a,b-a)
	else
		p =  M \ [imag(w);-real(w)] 
		cen = (a+b)/2 - 1im*p[1]*d1
		return Circle{Complex}(cen,abs(a-cen))
	end
end

point(C::Circle,t::Real) = C.center + C.radius*exp(2im*pi*t)
arclength(C::Circle) = 2π*C.radius

for T in [Complex,Polar,Spherical]
	struct Line{T} <: AbstractCurve 
		base::T 
		direction::T
		Line{T}(base,dir) = new(T(base),T(sign(dir)))
	end
end
Line(z::Complex,w::Complex) = Line{Complex}(z,w)
Line(z::Polar,w::Polar) = Line{Polar}(z,w)
Line(z::Spherical,w::Spherical) = Line{Spherical}(z,w)

point(L::Line,t::Real) = L.base + (2t-1)/(t-t^2)*L.direction
arclength(::Line) = Inf

# struct Line{Complex{T<:Real}} <: AbstractCurve 
# 	base::Complex{T}
# 	direction::Complex{T}
# 	Line{Complex{T}}(base,dir) = new(base,sign(dir)) 
# end

# struct Line{Polar{T<:Real}} <: AbstractCurve 
# 	base::Complex{T}
# 	direction::Complex{T}
# 	Line{Complex{T}}(base,dir) = new(base,sign(dir)) 
# end


struct Segment{T<:Number} <: AbstractCurve
	base::T 
	delta::T
	#Segment{T}(a,b) = new(T(a),T(b-a))
end

Segment(a::Real,b::Real) = Segment(Complex(a,0),Complex(b,0))
Segment(a::Number,b::Number) = Segment(promote(a,b)...)
Segment(a::T,b::T) where T = Segment{T}(a,b)
Segment{T}(base::T;delta::T=zero(T)) where T<:Number = Segment{T}(base,base+delta)

arclength(c::Segment) = abs(c.delta)
point(c::Segment,t::Real) = c.base + t*c.delta 

+(c::Segment,z::Number) = Segment()

show(io::IO,c::Segment) = "Segment from $(c.base) to $(c.base+c.delta)"

