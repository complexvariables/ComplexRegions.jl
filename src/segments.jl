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

# Complex type converters
for ctype in [:Spherical,:Polar,:Complex]
	@eval begin 
		function $ctype(S::Segment{T}) where T<:AnyComplex 
			Segment($ctype(S.za),$ctype(S.zb))
		end	
	end
end

# Required methods
arclength(S::Segment) = abs(S.zb-S.za)
point(S::Segment,t::Real) = (1-t)*S.za + t*S.zb
(C::Segment)(t::Real) = point(C,t)
arg(S::Segment,z::Number) = (real(z) - real(S.za)) / (real(S.zb) - real(S.za))
tangent(S::Segment,t::Real) = tangent(S)
tangent(S::Segment) = sign(S.zb-S.za)

# Other methods
isfinite(::Segment) = true
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
sign(S::Segment) = tangent(S)

function isapprox(S1::Segment,S2::Segment;tol=DEFAULT[:tol])
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

plotdata(S::Segment{T}) where T<:Union{Complex,Polar} = [S.za,S.zb]
