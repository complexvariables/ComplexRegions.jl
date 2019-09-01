# Type  
"""
	(type) Segment{T<:AnyComplex} in the complex plane 

Each `Segment` type is parameterized according to the common type of its complex input arguments. 
"""
struct Segment{T<:AnyComplex} <: AbstractCurve 
	za::T 
	zb::T 
	function Segment{T}(a,b) where T<:AnyComplex
		@assert isfinite(a) && isfinite(b)
		new(a,b)
	end
end

# Untyped constructor
"""
	Segment(a,b)

Consruct a segment that starts at value `a` and ends at `b`. 
"""
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
""" 
	arg(S::Segment,z) 

Find the parameter argument `t` such that `S(t)==z` is true. 

This gives undefined results if `z` is not actually on the segment. 
"""
arg(S::Segment,z::Number) = (real(z) - real(S.za)) / (real(S.zb) - real(S.za))
tangent(S::Segment,t::Real) = S.zb - S.za
unittangent(S::Segment,t::Real=0) = sign(S.zb-S.za)

+(S::Segment,z::Number) = Segment(S.za+z,S.zb+z)
-(S::Segment) = Segment(-S.za,-S.zb)
*(S::Segment,z::Number) = Segment(S.za*z,S.zb*z)

"""
	inv(S) 
Invert the segment `S` through the origin. In general the inverse is an `Arc`, though the result is a `Segment` if `S` would pass through the origin when extended.
"""
function inv(S::Segment) 
	w = 1 ./ point(S,[0,0.5,1])
	Arc(w...)
end

# Other methods
isfinite(::Segment) = true
conj(S::Segment) = Segment(conj(S.za),conj(S.zb))
reverse(S::Segment) = Segment(S.zb,S.za)

sign(S::Segment) = unittangent(S)

"""
	isapprox(S1::Segment,S2::Segment; tol=<default>) 
	S1 ≈ S2 

Determine if `S1` and `S2` represent the same segment, irrespective of the type or values of its parameters. Identity is determined by agreement within `tol`, which is interpreted as the weaker of absolute and relative differences.
"""
function isapprox(S1::Segment,S2::Segment;tol=DEFAULT[:tol])
	return isapprox(S1.za,S2.za,rtol=tol,atol=tol) &&
		isapprox(S1.zb,S2.zb,rtol=tol,atol=tol)
end

""" 
	isleft(z,S::Segment) 

Determine whether the number `z` lies "to the left" of segment `S`. This means that the angle it makes with `tangent(S)` is in the interval (0,π).

Note that `isleft` and `isright` are *not* logical opposites; a point on the (extended) segment should give `false` in both cases.
"""
function isleft(z::Number,S::Segment) 
	a,b = S.za,S.zb
	(real(b)-real(a)) * (imag(z)-imag(a)) > (real(z)-real(a)) * (imag(b)-imag(a))
end
""" 
	isright(z,S::Segment) 

Determine whether the number `z` lies "to the right" of segment `S`. This means that the angle it makes with `tangent(S)` is in the interval (-π,0).

Note that `isleft` and `isright` are *not* logical opposites; a point on the (extended) segment should give `false` in both cases.
"""
function isright(z::Number,S::Segment) 
	a,b = S.za,S.zb
	(real(b)-real(a)) * (imag(z)-imag(a)) < (real(z)-real(a)) * (imag(b)-imag(a))
end

""" 
	dist(z,S::Segment) 

Compute the distance from number `z` to the segment `S`. 
"""
dist(z::Number,S::Segment) = abs(z - closest(z,S))
""" 
	closest(z,S::Segment) 

Find the point on segment `S` that lies closest to `z`.
"""
function closest(z::Number,S::Segment) 
	# translate and rotate segment to positive Re axis
	d = S.zb-S.za
	s = sign(d)
	ζ = (z-S.za)/s
	S.za + s*min( max(real(ζ),0), abs(d) ) 
end
""" 
	reflect(z,S::Segment) 

Reflect the value `z` across the extension of segment `S` to a line. (For reflection of a segment through a point, use translation and negation.)
"""
reflect(z::Number,S::Segment) = reflect(z,Line(S.za,S.zb))

# Display methods
function show(io::IO,S::Segment{T}) where {T}
	print(IOContext(io,:compact=>true),"Segment(",point(S,0),",",point(S,1),")")
end

function show(io::IO,::MIME"text/plain",S::Segment{T}) where {T}
	print(io,"Segment{$T} in the complex plane:\n   from (",point(S,0),") to (",point(S,1),")")
end

plotdata(S::Segment{T}) where T<:Union{Complex,Polar} = [S.za,S.zb]
