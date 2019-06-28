abstract type AbstractCurve end

point(c::AbstractCurve,t::Real) = @error "No point() method defined for type $(typeof(c))"
point(c::AbstractCurve,t::AbstractArray{T}) where T<:Real = [point(c,t) for t in t]
start(c::AbstractCurve) = point(c,0.0)
stop(c::AbstractCurve) = point(c,1.0)
arclength(c::AbstractCurve) = @error "No arclength() method defined for type $(typeof(c))"

abstract type AbstractClosedCurve <: AbstractCurve end

#
# Circle
#

# Type and constructors
struct Circle{T<:AnyComplex} <: AbstractClosedCurve 
	center::T
	radius::Float64
end
Circle(z::AnyComplex,r::Real) = Circle{typeof(z)}(z,r)
Circle(z::Number,r::Real) = Circle(complex(float(z)),r)

# Construction by three points
function Circle(a::Number,b::Number,c::Number) 
	a,b,c = promote( complex.(float.([a,b,c]))... )
	Circle(a,b,c)
end
function Circle(a::T,b::T,c::T) where {T<:AnyComplex}
	isinf(a) && return Line(b,c)
	isinf(b) && return Line(a,c)
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
		return Circle{T}(cen,abs(a-cen))
	end
end

# Required methods
point(C::Circle,t::Real) = C.center + C.radius*exp(2im*pi*t)
arclength(C::Circle) = 2π*C.radius
(C::Circle)(t::Real) = point(C,t)

# Other methods
+(C::Circle,z::Number) = Circle(C.center+z,C.radius)
+(z::Number,C::Circle) = Circle(C.center+z,C.radius)
-(C::Circle) = Circle(-C.center,C.radius)
-(C::Circle,z::Number) = Circle(C.center-z,C.radius)
-(z::Number,C::Circle) = Circle(z-C.center,C.radius)
*(C::Circle,z::Number) = Circle(C.center*z,C.radius*abs(z))
*(z::Number,C::Circle) = Circle(C.center*z,C.radius*abs(z))
/(C::Circle,z::Number) = Circle(C.center/z,C.radius/abs(z))

function show(io::IO,C::Circle)
	print(IOContext(io,:compact=>true),"Circle(",C.center,",",C.radius,")")
end
function show(io::IO,::MIME"text/plain",C::Circle{T}) where {T}
	print(io,"Circle{$T} in the complex plane:\n   centered at (",C.center,") with radius ",C.radius)
end

#
# Line
# 

# Type and constructors
# Default constructor is for the line through two points
struct Line{T<:AnyComplex} <: AbstractClosedCurve 
	base::T 
	direction::T
	Line{T}(a,b) where T = new(a,b-a)
end
function Line(a::Number,b::Number)
	a,b = promote(complex(float(a)),complex(float(b)))
	Line{typeof(a)}(a,b)
end

# Construct by giving a direction keyword
Line(z::Number;direction) = Line(z,z+direction)

# Required methods
arclength(::Line) = Inf
function point(L::Line,t::Real)
	if t==0
		z = Polar(Inf,angle(-L.direction))
	elseif t==1
		z = Polar(Inf,angle(L.direction))
	else
		z = L.base + (2t-1)/(t-t^2)*L.direction
	end
	return z 
end
(C::Line)(t::Real) = point(C,t)

# Other methods
+(L::Line,z::Number) = Line(L.base+z,direction=L.direction)
+(z::Number,L::Line) = Line(L.base+z,direction=L.direction)
-(L::Line) = Line(-L.base,direction=-L.direction)
-(L::Line,z::Number) = Line(L.base-z,direction=L.direction)
-(z::Number,L::Line) = Line(z-L.base,direction=-L.direction)
*(L::Line,z::Number) = Line(L.base*z,direction=L.direction*z)
*(z::Number,L::Line) = Line(L.base*z,direction=L.direction*z)
/(L::Line,z::Number) = Line(L.base/z,direction=L.direction/z)
slope(L::Line) = tan(angle(L.direction))

function show(io::IO,L::Line)
	print(IOContext(io,:compact=>true),"Line(...",L(0.5),"...",L((sqrt(5)-1)/2),"...)")
end
function show(io::IO,::MIME"text/plain",L::Line{T}) where {T}
	print(io,"Line{$T} in the complex plane:\n   through (",L.base,") parallel to (",L.direction,")")
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
		delta = mod(angle(β/α)/(2π),1)
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
		delta = mod(angle(β/α)/(2π),1)
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

function show(io::IO,A::Arc{T}) where {T}
	print(IOContext(io,:compact=>true),"Arc(",A(0.0),"...",A(1.0),")")
end
function show(io::IO,::MIME"text/plain",A::Arc{T}) where {T}
	print(io,"Arc{$T} in the complex plane:\n   fraction ",A.delta," of (",A.circle,") starting at ",A.start)
end


# 
# Segment 
# 

# Type  
struct Segment{T<:AnyComplex} <: AbstractCurve 
	line::Line{T} 
	start::Float64  # specified as Line parameter value
	stop::Float64 
end

# Untyped constructor
function Segment(L::Line{T},start::Real,stop::Real) where T<:AnyComplex
	Segment{T}(L,Float64(start),Float64(stop))
end

# Construct from 2 points
function Segment(a::Number,b::Number) 
	a,b = promote(complex(float(a)),b)
	L = Line(a,b)
	tf = (sqrt(5)-1)/2  # based on Line implementation
	Segment(L,0.5,tf)
end

# Required methods
arclength(S::Segment) = abs(point(S,1) - point(S,0))
point(S::Segment,t::Real) = (1-t)*S.line(S.start) + t*S.line(S.stop)
(C::Segment)(t::Real) = point(C,t)

# Other methods
+(S::Segment,z::Number) = Segment(+(S.line,z),S.start,S.stop)
+(z::Number,S::Segment) = Segment(+(z,S.line),S.start,S.stop)
-(S::Segment,z::Number) = Segment(-(S.line,z),S.start,S.stop)
-(z::Number,S::Segment) = Segment(-(z,S.line),S.start,S.stop)
*(S::Segment,z::Number) = Segment(*(S.line,z),S.start,S.stop)
*(z::Number,S::Segment) = Segment(*(z,S.line),S.start,S.stop)
/(S::Segment,z::Number) = Segment(/(S.line,z),S.start,S.stop)
-(S::Segment) = Segment(-S.line,S.start,S.stop)

function show(io::IO,S::Segment{T}) where {T}
	print(IOContext(io,:compact=>true),"Segment(",point(S,0),",",point(S,1),")")
end
function show(io::IO,::MIME"text/plain",S::Segment{T}) where {T}
	print(io,"Segment{$T} in the complex plane:\n   from (",point(S,0),") to (",point(S,1),")")
end
