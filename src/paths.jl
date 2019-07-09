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

point(c::AbstractPath,t::AbstractArray{T}) where T<:Real = [point(c,t) for t in t]

eltype(::Type{AbstractPath}) = AbstractCurve 
length(p::AbstractPath) = length(curve(p))
getindex(p::AbstractPath,k) = curve(p,k)
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

conj(p::AbstractPath) = typeof(p)(conj.(curve(p)))
reverse(p::AbstractPath) = typeof(p)(reverse(reverse.(curve(p))))
+(p::AbstractPath,z::Number) = typeof(p)([c+z for c in curve(p)])
+(z::Number,p::AbstractPath) = typeof(p)([z+c for c in curve(p)])
-(p::AbstractPath) = typeof(p)([-c for c in curve(p)])
-(p::AbstractPath,z::Number) = typeof(p)([c-z for c in curve(p)])
-(z::Number,p::AbstractPath) = typeof(p)([z-c for c in curve(p)])
*(p::AbstractPath,z::Number) = typeof(p)([c*z for c in curve(p)])
*(z::Number,p::AbstractPath) = typeof(p)([z*c for c in curve(p)])
/(p::AbstractPath,z::Number) = typeof(p)([c/z for c in curve(p)])
isbounded(p::AbstractPath) = all( isfinite.(vertex(p)) )

function show(io::IO,P::AbstractPath)
	print(IOContext(io,:compact=>true),typeof(P)," with ",length(P)," segments") 
end
function show(io::IO,::MIME"text/plain",P::AbstractPath) 
	print(io,typeof(P)," with ",length(P)," segments")
end

abstract type AbstractClosedPath <: AbstractPath end

curve(p::AbstractClosedPath,k::Integer) = curve(p)[mod(k-1,length(p))+1]
vertex(P::AbstractClosedPath,k::Integer) = point(curve(P,k),0)
vertex(P::AbstractClosedPath) = [ vertex(P,k) for k = 1:length(P) ]

#
# Concrete implementations
#

# Path
struct Path <: AbstractPath 
	curve
	arclen
	breakindex
	function Path(p::AbstractVector;tol::Real=1e-13)
		n = length(p)
		for k = 1:n-1
			@assert p[k] isa AbstractCurve
			@assert isapprox(point(p[k],1.0),point(p[k+1],0.0),rtol=tol,atol=tol) "Curve endpoints do not match for pieces $(k) and $(k+1)"
		end
		@assert p[end] isa AbstractCurve
		len = [arclength(c) for c in p]
		s = cumsum(len)
		#new(p,len,s[1:end-1]/s[end])
		new(p,len,(1:n-1)/n)
	end
end
Path(c::AbstractCurve) = Path([c])
promote_rule(Curve,Path) = Path

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
promote_rule(ClosedCurve,ClosedPath) = ClosedPath

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

function rectangle(xlim::AbstractVector,ylim::AbstractVector)  
	x = [xlim[1],xlim[2],xlim[2],xlim[1],xlim[1]]
	y = [ylim[1],ylim[1],ylim[2],ylim[2],ylim[1]]
	Polygon( [Segment(complex(x[k],y[k]),complex(x[k+1],y[k+1])) for k in 1:4] )
end
rectangle(z1::AnyComplex,z2::AnyComplex) = rectangle([real(z1),real(z2)],[imag(z1),imag(z2)])
rectangle(z1::Number,z2::Number) = rectangle(promote(float(z1),float(z2))...)