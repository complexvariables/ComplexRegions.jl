# Type and constructors
# Default constructor is for the line through two points
"""
(type) Line in the complex plane 
"""
struct Line{T<:AnyComplex} <: AbstractClosedCurve 
	base::T 
	direction::T
	Line{T}(a,b) where T = new(a,sign(b-a))
end
"""
	Line(a,b)
	Line(a,direction=z)

Construct a line, either through the two given points, or by giving one point and a direction value `z` whose complex sign is parallel to the line. 
"""
function Line(a::Number,b::Number)
	a,b = promote(complex(float(a)),complex(float(b)))
	Line{typeof(a)}(a,b)
end

# Construct by giving a keyword
Line(z::Number;direction) = Line(z,z+direction)
#Line(z::Number;angle) = Line(z,z+exp(complex(0,angle)))

# Complex type converters
for ctype in [:Spherical,:Polar,:Complex]
	@eval begin 
		function $ctype(L::Line{T}) where T<:AnyComplex 
			Line($ctype(L.base),$ctype(L.base+L.direction))
		end
	end
end	

# Required methods
arclength(::Line) = Inf
point(L::Line,t::Real) = L.base + (2t-1)/(t-t^2)*L.direction
(C::Line)(t::Real) = point(C,t)
function arg(L::Line,z::Number)
	isinf(z) && return float(0)
	del = z-L.base
	if abs(real(del)) > abs(imag(del)) 
		del = real(del) 
		α = real(L.direction)
	else
		del = imag(del) 
		α = imag(L.direction)
	end
	for t in realroots(del,2α-del,-α)
		0-eps(del) ≤ t ≤ 1+eps(del) && return t 
	end
	return []
end
tangent(L::Line) = L.direction
tangent(L::Line,t::Real) = tangent(L)

# Other methods
isbounded(::Line) = false
slope(L::Line) = imag(L.direction)/real(L.direction)
conj(L::Line) = Line(conj(L.base),direction=conj(L.direction))
reverse(L::Line) = Line(L.base,direction=-L.direction)
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
inv(L::Line) = 1/L

function isapprox(L1::Line,L2::Line;tol=1e-12)
	dz = L1.base - L2.base
	w1,w2 = L1.direction,L2.direction
	return isapprox(real(w1)*imag(w2),imag(w1)*real(w2),rtol=tol,atol=tol) &&
		isapprox(real(w1)*imag(dz),imag(w1)*real(dz),rtol=tol,atol=tol)
end

isleft(z::Number,L::Line) = angle((z-L.base)/L.direction) > 0 

dist(z::Number,L::Line) = imag( (z-L.base)/L.direction )
function closest(z::Number,L::Line) 
	s = L.direction
	L.base + real((z-L.base)/s)*s 
end	
function reflect(z::Number,L::Line) 
	ζ = z - L.base 
	L.base + L.direction*conj(ζ/L.direction)
end

function show(io::IO,L::Line)
	print(IOContext(io,:compact=>true),"Line(...",L(0.5),"...",L((sqrt(5)-1)/2),"...)")
end
function show(io::IO,::MIME"text/plain",L::Line{T}) where {T}
	print(io,"Line{$T} in the complex plane:\n   through (",L.base,") parallel to (",L.direction,")")
end

plotdata(L::Line{T}) where T<:Union{Complex,Polar} = point(L,[0.1,0.9])
