abstract type AbstractPath end

# Required methods
"""
	curves(P::AbstractPath)
Return an array of the curves that make up the path `P`. 
"""
curves(p::AbstractPath) = @error "No curves() method defined for type $(typeof(p))"

# Default implementations
"""
	curve(P::AbstractPath,k::Integer)
Return the `k`th curve in the path `P`. 
"""
curve(p::AbstractPath,k::Integer) = curves(p)[k]

"""
	vertex(P::AbstractPath,k::Integer) 
Return the `k`th vertex of the path `P`.
"""
function vertex(P::AbstractPath,k::Integer) 
	C = curves(P)
	n = length(C)
	if 1 ≤ k ≤ n
		point(C[k],0)
	elseif k==n+1
		point(C[n],1)
	else
		throw(BoundsError(P,k))
	end
end

"""
	vertices(P::AbstractPath)
Return an array of the vertices (endpoints of the curves) of the path `P`. The length is one greater than the number of curves in `P`.
"""
function vertices(P::AbstractPath) 
	[ vertex(P,k) for k = 1:length(P)+1]
end

""" 
	isfinite(P::AbstractPath) 
Return `true` if the path is bounded in the complex plane (i.e., does not pass through infinity).
"""
isfinite(p::AbstractPath) = all(isfinite(s) for s in curves(p))

# iteration interface
eltype(::Type{AbstractPath}) = AbstractCurve 
length(p::AbstractPath) = length(curves(p))
getindex(p::AbstractPath,k) = curve(p,k)
iterate(p::AbstractPath,state=1) = state > length(curves(p)) ? nothing : (p[state], state+1)

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

conj(p::AbstractPath) = typeof(p)(conj.(curves(p)))

reverse(p::AbstractPath) = typeof(p)(reverse(reverse.(curves(p))))

+(p::AbstractPath,z::Number) = typeof(p)([c+z for c in curves(p)])
+(z::Number,p::AbstractPath) = typeof(p)([z+c for c in curves(p)])

-(p::AbstractPath) = typeof(p)([-c for c in curves(p)])
-(p::AbstractPath,z::Number) = typeof(p)([c-z for c in curves(p)])
-(z::Number,p::AbstractPath) = typeof(p)([z-c for c in curves(p)])

*(p::AbstractPath,z::Number) = typeof(p)([c*z for c in curves(p)])
*(z::Number,p::AbstractPath) = typeof(p)([z*c for c in curves(p)])

/(p::AbstractPath,z::Number) = typeof(p)([c/z for c in curves(p)])
/(z::Number,p::AbstractPath) = z*inv(p)
inv(p::AbstractPath) = typeof(p)([inv(c) for c in curves(p)])

"""
	isapprox(P1::AbstractPath,R2::AbstractPath; tol=<default>)
	P1 ≈ P2       (type "\\approx" followed by tab key)
Determine whether `P1` and `P2` represent the same path, up to tolerance `tol`, irrespective of the parameterization of its curves.
"""
function isapprox(P1::AbstractPath,P2::AbstractPath;tol=DEFAULT[:tol]) 
	if length(P1) != length(P2) 
		return false
	else
		c1,c2 = curves(P1),curves(P2)
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
	k = argmin( [dist(z,s) for s in curves(P)] )
	closest(z,side(P,k))
end

intersect(P::AbstractPath,C::AbstractCurve) = ∪( [intersect(s,C) for s in curves(P)]...  )
intersect(P1::AbstractPath,P2::AbstractPath) = ∪( [intersect(P1,s) for s in curves(P2)]...  )

function show(io::IO,P::AbstractPath)
	print(IOContext(io,:compact=>true),typeof(P)," with ",length(P)," curves") 
end
function show(io::IO,::MIME"text/plain",P::AbstractPath) 
	print(io,typeof(P)," with ",length(P)," curves")
end

plotdata(P::AbstractPath) = vcat(plotdata.(P)...)

abstract type AbstractClosedPath <: AbstractPath end

"""
	curve(P::AbstractClosedPath,k::Integer)
Return the `k`th curve in the path `P`. The index is applied circularly; e.g, if the closed path has n curves, then ...,1-n,1,1+n,... all refer to the first curve. 
"""
curve(p::AbstractClosedPath,k::Integer) = curves(p)[mod(k-1,length(p))+1]

"""
	vertex(P::AbstractPath,k::Integer) 
Return the `k`th vertex of the path `P`. The index is applied circularly; e.g, if the closed path has n curves, then ...,1-n,1,1+n,... all refer to the first vertex. 
"""
vertex(P::AbstractClosedPath,k::Integer) = point(curve(P,k),0)

"""
	vertices(P::AbstractClosedPath)
Return an array of the unique vertices (endpoints of the curves) of the closed path `P`. The length is equal the number of curves in `P`, i.e., the first/last vertex is not duplicated.
"""
vertices(P::AbstractClosedPath) = [ vertex(P,k) for k = 1:length(P) ]

function sideargs(p::AbstractClosedPath,t) 
	n = length(p)
	return 1+mod(floor(Int,t),n), t%1
end

function winding(P::AbstractClosedPath,z::Number)
	# Integrate around the boundary
	w = 0
	for s in P
		f = t -> imag(tangent(s,t)/(point(s,t)-z))
		w += intadapt(f,0,1,1e-4)
	end
	return round(Int,w/(2π))
end

isinside(z::Number,P::AbstractClosedPath) = winding(P,z) != 0 
isoutside(z::Number,P::AbstractClosedPath) = winding(P,z) == 0 

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

curves(p::Path) = p.curve 
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
		zi,zf = point(q,0),point(q,length(q))
		@assert isapprox(zi,zf,rtol=tol,atol=tol) || (isinf(zi) && isinf(zf)) "Path endpoints do not match"
		new(p)
	end
end
"""
	ClosedPath(c::AbstractVector; tol=<default>)
	ClosedPath(P::Path; tol=<default>)
Given a vector `c` of curves, or an existing path, construct a closed path. The path is checked for continuity (to tolerance `tol`) at all of the vertices. 
"""
ClosedPath(c::AbstractCurve) = ClosedPath([c])
ClosedPath(c::AbstractClosedPath) = c
ClosedPath(p::Path;kw...) = ClosedPath(p.curve;kw...)

curves(p::ClosedPath) = p.curve 
arclength(p::ClosedPath) = sum(arclength(c) for c in p)
(p::ClosedPath)(t) = point(p,t)

# Find a circle that fully encloses all the finite vertices and some points of a path.
function enclosing_circle(p::AbstractPath,expansion=2)
	#z = [vertices(p);point(p,0:1/50:1)]
	z = vertices(p)
	return Circle(enclosing_circle(filter(isfinite,z),expansion)...)
end

include("polygons.jl")
