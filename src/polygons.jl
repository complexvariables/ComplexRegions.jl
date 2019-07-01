abstract type AbstractCircularPolygon <: AbstractClosedPath end

struct CircularPolygon <: AbstractCircularPolygon
	side
	arclen
	breakindex
	function CircularPolygon(p::AbstractVector,arclen,breakindex)
		# Assumes continuity and closure have been checked previously
		valid = isa.(p,Union{Arc,Segment})
		@assert all(valid) "All sides must be an Arc or a Segment"
		new(p,arclen,breakindex)
	end
end

CircularPolygon(p::ClosedPath) = CircularPolygon(p.curve,p.arclen,p.breakindex)
CircularPolygon(p::Path;kw...) = CircularPolygon(ClosedPath(p;kw...))
function CircularPolygon(p::AbstractVector{T};kw...) where T<:AbstractCurve 
	CircularPolygon(ClosedPath(p;kw...))
end

curve(p::CircularPolygon) = p.side 
breakindex(p::CircularPolygon) = p.breakindex
arclength(p::CircularPolygon) = sum(p.arclen)
(p::CircularPolygon)(t::Real) = point(p,t)

# function show(io::IO,L::Line)
# 	print(IOContext(io,:compact=>true),"Line(...",L(0.5),"...",L((sqrt(5)-1)/2),"...)")
# end
# function show(io::IO,::MIME"text/plain",L::Line{T}) where {T}
# 	print(io,"Line{$T} in the complex plane:\n   through (",L.base,") parallel to (",L.direction,")")
# end

# 
# Polygon 
# 
abstract type AbstractPolygon <: AbstractCircularPolygon end

# Type 
struct Polygon <: AbstractPolygon
	side
	arclen
	breakindex
	function Polygon(p::AbstractVector{T},arclen,breakindex) where T<:AbstractCurve
		# Assumes continuity and closure have been checked previously
		valid = isa.(p,Segment)
		@assert all(valid) "All sides must be an Arc or a Segment"
		new(p,arclen,breakindex)
	end
end

Polygon(p::ClosedPath) = Polygon(p.curve,p.arclen,p.breakindex)
Polygon(p::Path;kw...) = Polygon(ClosedPath(p;kw...))
function Polygon(p::AbstractVector{T};kw...) where T<:AbstractCurve 
	Polygon(ClosedPath(p;kw...))
end
function Polygon(v::AbstractVector{T};kw...) where T<:Number 
	n = length(v)
	p = [Segment(v[j],v[mod(j,n)]) for j = 1:n]
	Polygon(ClosedPath(p;kw...))
end

curve(p::Polygon) = p.side 
breakindex(p::Polygon) = p.breakindex
arclength(p::Polygon) = sum(p.arclen)
(p::Polygon)(t::Real) = point(p,t)

# function show(io::IO,S::Segment{T}) where {T}
# 	print(IOContext(io,:compact=>true),"Segment(",point(S,0),",",point(S,1),")")
# end
# function show(io::IO,::MIME"text/plain",S::Segment{T}) where {T}
# 	print(io,"Segment{$T} in the complex plane:\n   from (",point(S,0),") to (",point(S,1),")")
# end
