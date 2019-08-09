abstract type AbstractPath end

# Required methods
"""
	curve(P::AbstractPath)
Return an array of the curves that make up the path `P`. 

	curve(P::AbstractPath,k::Integer)
Return the `k`th curve in the path `P`. 
"""
curve(p::AbstractPath) = @error "No curve() method defined for type $(typeof(p))"

# Default implementations
curve(p::AbstractPath,k::Integer) = curve(p)[k]
"""
	vertex(P::AbstractPath)
Return an array of the vertices (endpoints of the curves) of the path `P`. The length is one greater than the number of curves in `P`.

	vertex(P::AbstractPath,k::Integer) 
Return the `k`th vertex of the path `P`.
"""
function vertex(P::AbstractPath,k::Integer) 
	C = curve(P)
	n = length(C)
	if 1 ≤ k ≤ n
		point(C[k],0)
	elseif k==n+1
		point(C[n],1)
	else
		throw(BoundsError(P,k))
	end
end
function vertex(P::AbstractPath) 
	[ vertex(P,k) for k = 1:length(P)+1]
end

""" 
	isfinite(P::AbstractPath) 
Return `true` if the path is bounded in the complex plane (i.e., does not pass through infinity).
"""
isfinite(p::AbstractPath) = all(isfinite(s) for s in curve(p))

# iteration interface
eltype(::Type{AbstractPath}) = AbstractCurve 
length(p::AbstractPath) = length(curve(p))
getindex(p::AbstractPath,k) = curve(p,k)
iterate(p::AbstractPath,state=1) = state > length(curve(p)) ? nothing : (p[state], state+1)

"""
	point(P::AbstractPath,t::Real)
	P(t)
Compute the point along path `P` at parameter value `t`. Values of `t` in [k,k+1] correspond to values in [0,1] along curve k of the path, for k = 1,2,...,length(P)-1. 

	point(P::AbstractPath,t::AbstractVector)
Vectorize the `point` method for path `P`. 
"""
point(p::AbstractPath,t::AbstractArray{T}) where T<:Real = [point(p,t) for t in t]

# This determines how to parse a parameter of a path. Overloaded later for the case of a closed path.
function sideargs(p::AbstractPath,t) 
	n = length(p)
	if (t < 0) || (t > n) 
		throw(BoundsError(p,t))
	end
	if t==n 
		return n,1
	else
		return 1+floor(Int,t),t%1 
	end
end

function point(p::AbstractPath,t::Real) 
	k,s = sideargs(p,t)
	point(curve(p,k),s)
end

"""
	tangent(P::AbstractPath,t::Real)
Compute the complex-valued tangent along path `P` at parameter value `t`. Values of `t` in [k,k+1] correspond to values in [0,1] along curve k of the path, for k = 1,2,...,length(P)-1. The result is not well-defined at an integer value of `t`. 
"""
function tangent(p::AbstractPath,t::Real)
	k,s = sideargs(p,t)
	tangent(curve(p,k),s)
end

"""
	unittangent(P::AbstractPath,t::Real)
Compute the complex-valued unit tangent along path `P` at parameter value `t`. Values of `t` in [k,k+1] correspond to values in [0,1] along curve k of the path, for k = 1,2,...,length(P)-1. The result is not well-defined at an integer value of `t`. 
"""
function unittangent(p::AbstractPath,t::Real)
	k,s = sideargs(p,t)
	unittangent(curve(p,k),s)
end 

"""
	normal(P::AbstractPath,t::Real)
Compute a complex-valued normal to path `P` at parameter value `t`. Values of `t` in [k,k+1] correspond to values in [0,1] along curve k of the path, for k = 1,2,...,length(P)-1. The result is not well-defined at an integer value of `t`. 
"""
function normal(p::AbstractPath,t::Real)
	k,s = sideargs(p,t)
	normal(curve(p,k),s)
end

"""
	conj(P::AbstractPath) 
Construct the complex conjugate of `P`. Note that this also reverses the orientation of a closed path. 
"""
conj(p::AbstractPath) = typeof(p)(conj.(curve(p)))

""" 
	reverse(P::AbstractPath) 
Construct a path identical to `P` except with opposite direction of parameterization.
"""
reverse(p::AbstractPath) = typeof(p)(reverse(reverse.(curve(p))))

"""
	P + z
	z + P 
Translate the path `P` by a number `z`. 
"""
+(p::AbstractPath,z::Number) = typeof(p)([c+z for c in curve(p)])
+(z::Number,p::AbstractPath) = typeof(p)([z+c for c in curve(p)])

"""
	P - z
Translate the path `P` by a number `-z`.

	-P 
	z - P 
Negate a path `P` (reflect through the origin), and optionally translate by a number `z`.
"""
-(p::AbstractPath) = typeof(p)([-c for c in curve(p)])
-(p::AbstractPath,z::Number) = typeof(p)([c-z for c in curve(p)])
-(z::Number,p::AbstractPath) = typeof(p)([z-c for c in curve(p)])

"""
	z*P 
	P*z 
Multiply the path `P` by real or complex number `z`; i.e., scale and rotate it about the origin.
"""
*(p::AbstractPath,z::Number) = typeof(p)([c*z for c in curve(p)])
*(z::Number,p::AbstractPath) = typeof(p)([z*c for c in curve(p)])

"""
	P/z 
Multiply the path `P` by the number `1/z`; i.e., scale and rotate it about the origin.

	z/P 
	inv(P) 
Invert the path `P` through the origin (and optionally multiply by the number `z`). 
"""
/(p::AbstractPath,z::Number) = typeof(p)([c/z for c in curve(p)])
/(z::Number,p::AbstractPath) = typeof(p)([z/c for c in curve(p)])
inv(p::AbstractPath) = typeof(p)([inv(c) for c in curve(p)])

"""
	isapprox(P1::AbstractPath,R2::AbstractPath; tol=<default>)
	P1 ≈ P2       (type "\\approx" followed by tab key)
Determine whether `P1` and `P2` represent the same path, up to tolerance `tol`, irrespective of the parameterization of its curves.
"""
function isapprox(P1,P2;tol=DEFAULT[:tol]) 
	if length(P1) != length(P2) 
		return false
	else
		c1,c2 = curve(P1),curve(P2)
		all( isapprox(c1[k],c2[k]) for k in eachindex(c1) )
	end
end
isapprox(::AbstractCurve,::AbstractPath;kw...) = false
isapprox(::AbstractPath,::AbstractCurve;kw...) = false

""" 
	dist(z,P::AbstractPath) 
Find the distance from the path `P` to the point `z`.
"""
dist(z::Number,P::AbstractPath) = minimum(dist(z,C) for C in P)

""" 
	closest(z,P::AbstractPath) 
Find the point on the path `P` that lies closest to `z`.
"""
function closest(z::Number,P::AbstractPath)
	k = argmin( [dist(z,s) for s in curve(P)] )
	closest(z,side(P,k))
end

intersect(P::AbstractPath,C::AbstractCurve) = ∪( [intersect(s,C) for s in curve(P)]...  )
intersect(P1::AbstractPath,P2::AbstractPath) = ∪( [intersect(P1,s) for s in curve(P2)]...  )

function show(io::IO,P::AbstractPath)
	print(IOContext(io,:compact=>true),typeof(P)," with ",length(P)," curves") 
end
function show(io::IO,::MIME"text/plain",P::AbstractPath) 
	print(io,typeof(P)," with ",length(P)," curves")
end

abstract type AbstractClosedPath <: AbstractPath end

"""
	curve(P::AbstractClosedPath,k::Integer)
Return the `k`th curve in the path `P`. The index is applied circularly; e.g, if the closed path has n curves, then ...,1-n,1,1+n,... all refer to the first curve. 
"""
curve(p::AbstractClosedPath,k::Integer) = curve(p)[mod(k-1,length(p))+1]
"""
	vertex(P::AbstractClosedPath)
Return an array of the unique vertices (endpoints of the curves) of the closed path `P`. The length is equal the number of curves in `P`, i.e., the first/last vertex is not duplicated.

	vertex(P::AbstractPath,k::Integer) 
Return the `k`th vertex of the path `P`. The index is applied circularly; e.g, if the closed path has n curves, then ...,1-n,1,1+n,... all refer to the first vertex. 
"""
vertex(P::AbstractClosedPath,k::Integer) = point(curve(P,k),0)
vertex(P::AbstractClosedPath) = [ vertex(P,k) for k = 1:length(P) ]

function sideargs(p::AbstractClosedPath,t) 
	n = length(p)
	return 1+mod(floor(Int,t),n), t%1
end

#
# Concrete implementations
#

# Path
"""
	(type) Path 
Generic implementation of an `AbstractPath`.
"""
struct Path <: AbstractPath 
	curve
	function Path(c::AbstractVector{T};tol::Real=DEFAULT[:tol]) where T<:AbstractCurve
		n = length(c)
		for k = 1:n-1
				@assert isapprox(point(c[k],1.0),point(c[k+1],0.0),rtol=tol,atol=tol) "Curve endpoints do not match for pieces $(k) and $(k+1)"
		end
		new(c)
	end
end
"""
	Path(c::AbstractVector; tol=<default>)
Given a vector `c` of curves, construct a path. The path is checked for continuity (to tolerance `tol`) at the interior vertices. 
"""
Path(c::AbstractCurve) = Path([c])

curve(p::Path) = p.curve 
"""
	arclength(P::Path)
Compute the arclength of the path `P`.
"""
arclength(p::Path) = sum(arclength(c) for c in p)
(p::Path)(t) = point(p,t)

# ClosedPath
"""
	(type) ClosedPath 
Generic implementation of an `AbstractClosedPath`.
"""
struct ClosedPath <: AbstractClosedPath 
	curve
	function ClosedPath(p::AbstractVector{T};tol=DEFAULT[:tol]) where T<:AbstractCurve
		q = Path(p)
		@assert isapprox(point(q,length(q)),point(q,0),rtol=tol,atol=tol) "Path endpoints do not match"
		new(p)
	end
end
"""
	ClosedPath(c::AbstractVector; tol=<default>)
	ClosedPath(P::Path; tol=<default>)
Given a vector `c` of curves, or an existing path, construct a closed path. The path is checked for continuity (to tolerance `tol`) at all of the vertices. 
"""
ClosedPath(c::AbstractCurve) = ClosedPath([c])
ClosedPath(p::Path;kw...) = ClosedPath(p.curve;kw...)

curve(p::ClosedPath) = p.curve 
arclength(p::ClosedPath) = sum(arclength(c) for c in p)
(p::ClosedPath)(t) = point(p,t)

# Still experimental; don't document.
function isleft(z::Number,P::AbstractClosedPath)
	# TODO: this isn't foolproof
	d = dist(z,P) 
	t = 0:d/10:1 
	if t[end]==1
		t = t[1:end-1]
	end
	isleft(z,Polygon(P.(t)))
end

# 
include("polygons.jl")

"""
	rectangle(xlim,ylim) 
Construct the rectangle defined by `xlim[1]`` < x < `xlim[2]`, `ylim[1]`` < y < `ylim[2]`, where z=complex(x,y).
"""
function rectangle(xlim::AbstractVector,ylim::AbstractVector)  
	x = [xlim[1],xlim[2],xlim[2],xlim[1],xlim[1]]
	y = [ylim[1],ylim[1],ylim[2],ylim[2],ylim[1]]
	Polygon( [Segment(complex(x[k],y[k]),complex(x[k+1],y[k+1])) for k in 1:4] )
end

""" 
	rectangle(z1,z2) 
Construct the rectangle whose opposing corners are the given complex values. 
"""
rectangle(z1::AnyComplex,z2::AnyComplex) = rectangle([real(z1),real(z2)],[imag(z1),imag(z2)])
rectangle(z1::Number,z2::Number) = rectangle(promote(float(z1),float(z2))...)

"""
	n_gon(n)
Construct a regular n-gon with vertices on the unit circle.
"""
function n_gon(n::Integer) 
	@assert n > 2 "Must have at least three vertices"
	Polygon( exp.(2im*pi*(0:n-1)/n) ) 
end