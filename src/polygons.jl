abstract type AbstractCircularPolygon <: AbstractClosedPath end

# Common methods
side = curve

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

# Other methods
function winding(z::Number,p::CircularPolygon)
	# Approximate by a thoughtfully chosen polygon. The problem is that while we ought to customize the discretized polygon by how close z is to the boundary, it's going to get slow to do that for each new point z.  
	#p = truncate(p)  # TODO: no infinities
	n = length(p)
	L = arclength(p)
	w = Vector{Any}(undef,n) 
	for (k,s) in enumerate(side(p))
		if s isa Segment
			w[k] = [s.za,s.zb]
		else
			m = min(10*n,ceil(Int,20*arclength(s)/L))
			w[k] = point(s,LinRange(0,1,m+1))
		end
	end
	return winding(z,Polygon(vcat(w...)))
end

# TODO truncate circular polygons


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

function truncate(p::Polygon) 
	# try to find a circle clear of the polygon
	v = filter(isfinite,vertex(p))
	length(v) == length(p) && return p   # nothing to do
	zc = sum(v)/length(v) 
	R = maximum(@. abs(v - zc))
	return truncate(p,Circle(zc,2*R))
end

function truncate(p::Polygon,c::Circle) 
	n = length(p)
	s,v = side(p),vertex(p)
	snew = Vector{Any}(undef,n)
	z_pre = NaN
	if isa(s[1],Ray) && isa(s[n],Ray)
		# recognize that first side is actually the return from an infinite vertex 
		z_pre = intersect(c,s[n])[1]
	end
	for k = 1:n	
		if !isa(s[k],Ray)
			snew[k] = s[k] 
		else
			# first of a pair? 
			if isnan(z_pre) 
				z_pre = intersect(c,s[k])[1]
				snew[k] = Segment(s[k].base,z_pre)
			else
				z_post = intersect(c,s[k])[1]
				snew[k] = [Arc(z_pre,z_post,center=c.center),Segment(z_post,s[k].base)]
				z_pre = NaN
			end
		end
	end
	return CircularPolygon(vcat(snew...))
end

# TODO: unpredictable results for points on the boundary
# TODO: maybe wrong for some unbounded polygons
# Ref Dan Sunday, http://geomalgorithms.com/a03-_inclusion.html

function winding(z::Number,p::Polygon)
	function crossing(y,s::Segment)
		if y ≥ imag(s(0))
			down = false
			up = y < imag(s(1))
		else
			up = false
			down = y ≥ imag(s(1))
		end
		return up,down 
	end

	# To do a lot of points, it's better to truncate before calling. But here we go. 
	if !isbounded(p) 
		return winding(z,truncate(p))
	end

	wind = 0
	x,y = real(z),imag(z)
	for s in p.side
		up,down = crossing(y,s)
		@debug (s,up,down,isleft(z,s))
		if up && isleft(z,s) 
			wind += 1
		elseif down && !isleft(z,s) 
			wind -= 1
		end
	end
	return wind
end

isleft(z::Number,p::Polygon) = winding(z,p) > 0
isright(z::Number,p::Polygon) = winding(z,p) < 0