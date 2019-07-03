abstract type AbstractCircularPolygon <: AbstractClosedPath end

# Common methods

function show(io::IO,P::AbstractCircularPolygon)
	print(IOContext(io,:compact=>true),typeof(P)," with ",length(P)," sides") 
end
function show(io::IO,::MIME"text/plain",P::AbstractCircularPolygon) 
	print(io,typeof(P)," with ",length(P)," sides")
end

#
# CircularPolygon
#

struct CircularPolygon <: AbstractCircularPolygon
	side
	arclen
	breakindex
	function CircularPolygon(p::AbstractVector,arclen,breakindex)
		# Assumes continuity and closure have been checked previously
		valid = isa.(p,Union{Arc,Segment,Ray})
		@assert all(valid) "All sides must be an Arc, Segment, or Ray"
		new(p,arclen,breakindex)
	end
end

# Constructors
CircularPolygon(p::ClosedPath) = CircularPolygon(p.curve,p.arclen,p.breakindex)
CircularPolygon(p::Path;kw...) = CircularPolygon(ClosedPath(p;kw...))
function CircularPolygon(p::AbstractVector{T};kw...) where T<:AbstractCurve 
	CircularPolygon(ClosedPath(p;kw...))
end

# Required methods
curve(p::CircularPolygon) = p.side 
breakindex(p::CircularPolygon) = p.breakindex
arclength(p::CircularPolygon) = sum(p.arclen)
(p::CircularPolygon)(t::Real) = point(p,t)

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
		valid = isa.(p,Union{Segment,Ray})
		@assert all(valid) "All sides must be a Segment or Ray"
		new(p,arclen,breakindex)
	end
end

# Constructors
Polygon(p::ClosedPath) = Polygon(p.curve,p.arclen,p.breakindex)
Polygon(p::Path;kw...) = Polygon(ClosedPath(p;kw...))

function Polygon(p::AbstractVector{T};kw...) where T<:AbstractCurve 
	Polygon(ClosedPath(p;kw...))
end

function Polygon(v::AbstractVector;kw...)
	n = length(v)
	p = Vector{Union{Segment,Ray}}(undef,n)
	for j = 1:n
		vthis = v[j]
		vnext = v[mod(j,n)+1]
		if isa(vthis,Tuple)
			if isa(vnext,Tuple)
				@error("Cannot have consecutive infinite vertices")
			else
				p[j] = Ray(vnext,vthis[2],true)
			end 
		else
			if isa(vnext,Tuple)
				p[j] = Ray(vthis,vnext[1])
			else
				p[j] = Segment(vthis,vnext)
			end
		end
	end
	@debug for c in p 
		@show c 
		@show (c(0),c(1))
	end
	return Polygon(ClosedPath(p;kw...))
end

# Required methods
curve(p::Polygon) = p.side 
breakindex(p::Polygon) = p.breakindex
arclength(p::Polygon) = sum(p.arclen)
(p::Polygon)(t::Real) = point(p,t)

# Display methods 
function show(io::IO,::MIME"text/plain",P::Polygon) 
	print(io,"Polygon with ",length(P)," vertices:")
	for v in vertex(P)
		print("\n   ")
		show(io,MIME("text/plain"),v)
	end
end

# Other methods
function angle(p::Polygon)
	# computes a turn angle in (-pi,pi]  (neg = left turn)
	turn(s1,s2) = π - mod2pi(angle(s2/s1)+π)
	s = sign.(p) 
	n = length(p) 
	v = vertex(p)
	θ = similar(real(s))
	for k = 1:n 
		θ[k] = π + turn(s[mod(k-2,n)+1],s[k])
		if isinf(v[k]) 
			θ[k] -= 2π
		end
	end
	# correct for possible clockwise orientation
	sum(θ) > 0 ? θ : -θ
end

# unpredictable results for points on the boundary
# wrong for unbounded polygons
# Ref Dan Sunday, http://geomalgorithms.com/a03-_inclusion.html

function winding(z::Number,p::Polygon)
	wind = 0
	v = vertex(p)
	x,y = real(z),imag(z)
	n = length(v)
	for j = 1:n 
		vnext = v[mod(j,n)+1] 
		if imag(v[j]) ≤ y
			# upward crossing and to the left of the side
			if imag(vnext) > y  && isleft(z,p.side[j])
				wind += 1
			end
		else 
			# downward crossing and to the right of the side
			if imag(vnext) ≤ y  && !isleft(z,p.side[j])
				wind -= 1
			end
		end
	end
	return wind
end
