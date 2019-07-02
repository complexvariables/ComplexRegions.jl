# Line

# Type and constructors
# Default constructor is for the line through two points
struct Line{T<:AnyComplex} <: AbstractClosedCurve 
	base::T 
	direction::T
	Line{T}(a,b) where T = new(a,sign(b-a))
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
slope(L::Line) = imag(L.direction)/real(L.direction)

+(L::Line,z::Number) = Line(L.base+z,direction=L.direction)
+(z::Number,L::Line) = Line(L.base+z,direction=L.direction)
-(L::Line) = Line(-L.base,direction=-L.direction)
-(L::Line,z::Number) = Line(L.base-z,direction=L.direction)
-(z::Number,L::Line) = Line(z-L.base,direction=-L.direction)
*(L::Line,z::Number) = Line(L.base*z,direction=L.direction*sign(z))
*(z::Number,L::Line) = Line(L.base*z,direction=L.direction*sign(z))
/(L::Line,z::Number) = Line(L.base/z,direction=L.direction/sign(z))
function /(z::Number,L::Line) 
	w = z./point(L,[0.25,0.5,0.75])
	Circle(w...)
end

function isapprox(L1::Line,L2::Line;tol=1e-12)
	dz = L1.base - L2.base
	w1,w2 = L1.direction,L2.direction
	return isapprox(real(w1)*imag(w2),imag(w1)*real(w2),rtol=tol,atol=tol) &&
		isapprox(real(w1)*imag(dz),imag(w1)*real(dz),rtol=tol,atol=tol)
end

dist(z::Number,L::Line) = imag( (z-L.base)/L.direction )
function closest(z::Number,L::Line) 
	s = L.direction
	L.base + real((z-L.base)/s)*s 
end	

function show(io::IO,L::Line)
	print(IOContext(io,:compact=>true),"Line(...",L(0.5),"...",L((sqrt(5)-1)/2),"...)")
end
function show(io::IO,::MIME"text/plain",L::Line{T}) where {T}
	print(io,"Line{$T} in the complex plane:\n   through (",L.base,") parallel to (",L.direction,")")
end

# 
# Segment 
# 

# Type  
struct Segment{T<:AnyComplex} <: AbstractCurve 
	base::T 
	delta::T 
	reverse::Bool
	function Segment{T}(a,b) where T<:Complex
		if isinf(a) || isinf(b)
			@error("Must use Polar or Spherical value to designate infinity")
		else
		 new(a,b-a,false)
		end
	end
	function Segment{T}(a,b) where T<:Union{Polar,Spherical}
		if isinf(a)
			if isinf(b)
				@error("Segment with two infinite endpoints is not defined")
			else
				new(b,a-b,true)
			end
		else
		 	new(a,b-a,false)
		end
	end
end

# Untyped constructor
function Segment(a::Number,b::Number) 
	a,b = promote(complex(float(a)),b)
	Segment{typeof(a)}(a,b)
end

# Required methods
arclength(S::Segment) = abs(S.delta)
function point(S::Segment,t::Real)
	# have to take care of infinite endpoint
	a,d = S.base,S.delta
	if S.reverse 
		t = 1-t
	end
	if isinf(d)
		a + t/(1-t)*sign(d) 
	else
		a + t*d
	end
end
(C::Segment)(t::Real) = point(C,t)

# Other methods
+(S::Segment,z::Number) = Segment(S(0)+z,S(1)+z)
+(z::Number,S::Segment) = Segment(S(0)+z,S(1)+z)
-(S::Segment,z::Number) = Segment(S(0)-z,S(1)-z)
-(z::Number,S::Segment) = Segment(z-S(0),z-S(1))
-(S::Segment) = Segment(-S(0),-S(1))
# these need to recompute the final parameter values
*(S::Segment,z::Number) = Segment(z*S(0),z*S(1))
*(z::Number,S::Segment) = Segment(z*S(0),z*S(1))
/(S::Segment,z::Number) = Segment(S(0)/z,S(1)/z)
function /(z::Number,S::Segment) 
	w = z./point(S,[0,0.5,1])
	Arc(w...)
end

function isapprox(S1::Segment,S2::Segment;tol=1e-12)
	return isapprox(S1(0.0),S2(0.0),rtol=tol,atol=tol) &&
		isapprox(S1(1.0),S2(1.0),rtol=tol,atol=tol)
end

dist(z::Number,S::Segment) = abs(z - closest(z,S))
function closest(z::Number,S::Segment) 
	# translate and rotate segment to positive Re axis
	a,d = S.base,S.delta
	s = sign(d)
	ζ = (z-a)/s
	if real(ζ) < 0
		a 
	elseif real(ζ) > abs(d)
		a+d 
	else 
		a + real(ζ)*s
	end 
end

function show(io::IO,S::Segment{T}) where {T}
	print(IOContext(io,:compact=>true),"Segment(",point(S,0),",",point(S,1),")")
end
function show(io::IO,::MIME"text/plain",S::Segment{T}) where {T}
	print(io,"Segment{$T} in the complex plane:\n   from (",point(S,0),") to (",point(S,1),")")
end
