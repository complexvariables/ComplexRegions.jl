abstract type AbstractPath end

# Required methods
curve(p::AbstractPath) = @error "No curve() method defined for type $(typeof(p))"
breakindex(p::AbstractPath) = @error "No breakindex() method defined for type $(typeof(p))"

# Methods in common
curve(p::AbstractPath,k::Integer) = curve(p)[k]
breakindex(p::AbstractPath,k::Integer) = breakindex(p)[k]
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

eltype(::Type{AbstractPath}) = AbstractCurve 
length(p::AbstractPath) = length(curve(p))
getindex(p::AbstractPath,k) = getindex(curve(p),k)
iterate(p::AbstractPath,state=1) = state > length(curve(p)) ? nothing : (p[state], state+1)

function point(p::AbstractPath,t::Real)
	@assert (0 ≤ t ≤ 1) "Parameter is out of the range [0,1]."
	offset = [-eps(); breakindex(p); 1]
	c = curve(p)
	j = findlast(t .> offset)
	offset[1] = 0;
	s = scalefrom(offset[j],offset[j+1],t) 
	point(c[j],s)
end

+(p::AbstractPath,z::Number) = typeof(p)([c+z for c in curve(p)])
+(z::Number,p::AbstractPath) = typeof(p)([z+c for c in curve(p)])
-(p::AbstractPath) = typeof(p)([-c for c in curve(p)])
-(p::AbstractPath,z::Number) = typeof(p)([c-z for c in curve(p)])
-(z::Number,p::AbstractPath) = typeof(p)([z-c for c in curve(p)])
*(p::AbstractPath,z::Number) = typeof(p)([c*z for c in curve(p)])
*(z::Number,p::AbstractPath) = typeof(p)([z*c for c in curve(p)])
/(p::AbstractPath,z::Number) = typeof(p)([c/z for c in curve(p)])

abstract type AbstractClosedPath <: AbstractPath end

function vertex(P::AbstractClosedPath,k::Integer) 
	C = curve(P)
	n = length(C)
	point(C[1+mod(k-1,n)],0)
end
function vertex(P::AbstractClosedPath) 
	[ vertex(P,k) for k = 1:length(P)]
end

#
# Concrete implementations
#

# Path
struct Path <: AbstractPath 
	curve
	arclen
	breakindex
	function Path(p::AbstractVector;tol::Real=1e-13)
		for k = 1:length(p)-1
			@assert p[k] isa AbstractCurve
			@assert isapprox(point(p[k],1.0),point(p[k+1],0.0),rtol=tol,atol=tol) "Curve endpoints do not match for pieces $(k) and $(k+1)"
		end
		@assert p[end] isa AbstractCurve
		len = [arclength(c) for c in p]
		s = cumsum(len)
		new(p,len,s[1:end-1]/s[end])
	end
end
Path(c::AbstractCurve) = Path([c])

curve(p::Path) = p.curve 
breakindex(p::Path) = p.breakindex
arclength(p::Path) = sum(p.arclen)
(p::Path)(t::Real) = point(p,t)

# ClosedPath
struct ClosedPath <: AbstractClosedPath 
	curve
	arclen
	breakindex
	function ClosedPath(p::AbstractVector;tol::Real=1e-13)
		q = Path(p)
		@assert isapprox(point(q,1.0),point(q,0.0),rtol=tol,atol=tol) "Path endpoints do not match"
		new(p,q.arclen,q.breakindex)
	end
end
ClosedPath(c::AbstractCurve) = ClosedPath([c])
ClosedPath(p::Path;kw...) = ClosedPath(p.curve;kw...)

curve(p::ClosedPath) = p.curve 
breakindex(p::ClosedPath) = p.breakindex
arclength(p::ClosedPath) = sum(p.arclen)
(p::ClosedPath)(t::Real) = point(p,t)