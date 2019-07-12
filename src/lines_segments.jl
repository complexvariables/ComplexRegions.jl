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
point(L::Line,t::Real) = L.base + (2t-1)/(t-t^2)*L.direction
(C::Line)(t::Real) = point(C,t)

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

plotdata(L::Line) = point(L,[0.1,0.9])

# 
# Segment 
# 

# Type  
struct Segment{T<:AnyComplex} <: AbstractCurve 
	za::T 
	zb::T 
	function Segment{T}(a,b) where T<:AnyComplex
		@assert isfinite(a) && isfinite(b)
		new(a,b)
	end
end

# Untyped constructor
function Segment(a::Number,b::Number) 
	a,b = promote(complex(float(a)),b)
	Segment{typeof(a)}(a,b)
end

# Required methods
arclength(S::Segment) = abs(S.zb-S.za)
point(S::Segment,t::Real) = (1-t)*S.za + t*S.zb
(C::Segment)(t::Real) = point(C,t)

# Other methods
isbounded(::Segment) = true
conj(S::Segment) = Segment(conj(S.za),conj(S.zb))
reverse(S::Segment) = Segment(S.zb,S.za)
+(S::Segment,z::Number) = Segment(S.za+z,S.zb+z)
+(z::Number,S::Segment) = Segment(S.za+z,S.zb+z)
-(S::Segment,z::Number) = Segment(S.za-z,S.zb-z)
-(z::Number,S::Segment) = Segment(z-S.za,z-S.zb)
-(S::Segment) = Segment(-S.za,-S.zb)
# these need to recompute the final parameter values
*(S::Segment,z::Number) = Segment(S.za*z,S.zb*z)
*(z::Number,S::Segment) = Segment(S.za*z,S.zb*z)
/(S::Segment,z::Number) = *(S,1/z)
function /(z::Number,S::Segment) 
	w = z./point(S,[0,0.5,1])
	Arc(w...)
end
inv(S::Segment) = 1/S
sign(S::Segment) = sign(S.zb-S.za)

function isapprox(S1::Segment,S2::Segment;tol=1e-12)
	return isapprox(S1.za,S2.za,rtol=tol,atol=tol) &&
		isapprox(S1.zb,S2.zb,rtol=tol,atol=tol)
end

dist(z::Number,S::Segment) = abs(z - closest(z,S))
function closest(z::Number,S::Segment) 
	# translate and rotate segment to positive Re axis
	d = S.zb-S.za
	s = sign(d)
	ζ = (z-S.za)/s
	S.za + s*min( max(real(ζ),0), abs(d) ) 
end
reflect(z::Number,S::Segment) = reflect(z,Line(S.za,S.zb))

function isleft(z::Number,S::Segment) 
	a,b = S.za,S.zb
	(real(b)-real(a)) * (imag(z)-imag(a)) > (real(z)-real(a)) * (imag(b)-imag(a))
end

# Display methods
function show(io::IO,S::Segment{T}) where {T}
	print(IOContext(io,:compact=>true),"Segment(",point(S,0),",",point(S,1),")")
end

function show(io::IO,::MIME"text/plain",S::Segment{T}) where {T}
	print(io,"Segment{$T} in the complex plane:\n   from (",point(S,0),") to (",point(S,1),")")
end

plotdata(S::Segment) = [S.za,S.zb]

#
# Ray
# 

# Type  
struct Ray{T<:AnyComplex} <: AbstractCurve 
	base::T 
	angle::AbstractFloat  
	reverse::Bool
	function Ray{T}(a,d,rev=false) where T<:AnyComplex
		new(a,mod2pi(d),rev)
	end
end

# Untyped constructor
function Ray(a::Number,d::Number,rev=false) 
	a = complex(float(a))
	Ray{typeof(a)}(a,float(d),rev)
end

# Required methods
arclength(R::Ray) = Inf
function point(R::Ray{T},t::Real) where T
	if R.reverse 
		t = 1-t 
	end
	# avoid NaNs 
	if t==0 
		R.base 
	elseif t==1 
		T(Inf)
	else
		R.base + t/(1-t)*exp(complex(0,R.angle))
	end	
end
(C::Ray)(t::Real) = point(C,t)

# Other methods
isbounded(::Ray) = false
conj(R::Ray) = Ray(conj(R.base),-R.angle,R.reverse)
reverse(R::Ray) = Ray(R.base,R.angle,!R.reverse)
+(R::Ray,z::Number) = Ray(R.base+z,R.angle,R.reverse)
+(z::Number,R::Ray) = Ray(R.base+z,R.angle,R.reverse)
-(R::Ray,z::Number) = Ray(R.base-z,R.angle,R.reverse)
-(z::Number,R::Ray) = Ray(z-R.base,R.angle,R.reverse)
-(R::Ray) = Ray(-R.base,mod2pi(R.angle+pi),R.reverse)
# these need to recompute the final parameter values
*(R::Ray,z::Number) = Ray(z*R.base,mod2pi(R.angle+sign(z)),R.reverse)
*(z::Number,R::Ray) = *(R,z)
/(R::Ray,z::Number) = *(R,1/z)
function /(z::Number,R::Ray) 
	w = z./point(R,[0,0.5,1])
	Arc(w...)
end
inv(R::Ray) = 1/R

function isapprox(R1::Ray,R2::Ray;tol=1e-12)
	return isapprox(R1.base,R2.base,rtol=tol,atol=tol) && (abs(mod2pi(R1.angle-R2.angle)) < tol)
end

dist(z::Number,R::Ray) = abs(z - closest(z,R))
function closest(z::Number,R::Ray) 
	# translate and rotate to positive Re axis
	s = exp(complex(0,R.angle))
	ζ = (z-R.base)/s
	R.base + max(real(ζ),0)*s
end

sign(R::Ray) = R.reverse ? -exp(complex(0,R.angle)) : exp(complex(0,R.angle))

function isleft(z::Number,R::Ray) 
	a,b = point(R,[0.2,0.8])  # accounts for reversal
	(real(b)-real(a)) * (imag(z)-imag(a)) > (real(z)-real(a)) * (imag(b)-imag(a))
end

# Display methods
function show(io::IO,R::Ray{T}) where {T}
	print(IOContext(io,:compact=>true),"Ray(",R.base,",",R.angle,",",R.reverse,")")
end

function show(io::IO,::MIME"text/plain",R::Ray{T}) where {T}
	if R.reverse 
		print(io,"Ray to ",R.base," at angle ",R.angle)
	else
		print(io,"Ray from ",R.base," at angle ",R.angle)
	end
end

plotdata(R::Line) = R.reverse ? [R(0.3),R.base] : [R.base,R(0.7)]
