abstract type AbstractCircularPolygon <: AbstractClosedPath end
abstract type AbstractPolygon <: AbstractCircularPolygon end

function show(io::IO,P::AbstractCircularPolygon)
	print(IOContext(io,:compact=>true),typeof(P)," with ",length(P)," sides")
end
function show(io::IO,::MIME"text/plain",P::AbstractCircularPolygon)
	print(io,typeof(P)," with ",length(P)," sides")
end

# Other methods
sides(p::AbstractCircularPolygon) = curves(p)
side(p::AbstractCircularPolygon,args...) = curve(p,args...)

function reverse(p::AbstractCircularPolygon)
	@assert sum(isinf.(vertices(p))) < 2 "Reversal is not a propoer polygon"
	typeof(p)(reverse(reverse.(sides(p))))
end

function winding(p::AbstractCircularPolygon,z::Number)
	if isinf(z)
		!isfinite(p) && @warn "Winding number is ill-defined for a boundary point"
		return 0
	end

	# truncate an unbounded path
	if !isfinite(p)
		C = enclosing_circle(p)
		while !isinside(z,C)
			C = Circle(C.center,2*C.radius)
		end
		return winding(truncate(p,C),z)
	end

	w = 0
	v = vertex(p,1)
	for s in p
		vnew = point(s,1)
		if s isa Segment
			w += angle((vnew-z)/(v-z))
		else  #if s isa Arc
			# move the branch cut to avoid the arc
			if abs(z-s.circle.center) < s.circle.radius
				# can use ray from the center to a point not on s
				u = point(s.circle,s.start-(1-s.delta)/2) # not on s
				w += angle((vnew-z)/(z-u)) - angle((v-z)/(z-u))
			else
				# can put center on the positive real axis
				u = s.circle.center
				w += angle((vnew-z)/(u-z)) - angle((v-z)/(u-z))
			end
		end
		v = vnew
	end
	if isnan(w)
		@warn "Point may be on the boundary"
		return 1
	end
	return round(Int,w/(2π))
end

#
# CircularPolygon
#
"""
	(type) CircularPolygon
Type for closed paths consisting entirely of arcs, segments, and rays.
"""
struct CircularPolygon <: AbstractCircularPolygon
	path
	function CircularPolygon(p::AbstractClosedPath)
		# Continuity and closure have been checked to make a closed path
		valid = isa.(curves(p),Union{Arc,Segment,Ray})
		@assert all(valid) "All sides must be an Arc, Segment, or Ray"
		new(p)
	end
end

# Constructors
"""
	CircularPolygon(p::AbstractPath; tol=<default>)
	CircularPolygon(p::AbstractVector; tol=<default>)
Construct a circular polygon from a (possibly closed) path, or from a vector of curves. The `tol` parameter is a tolerance used when checking continuity and closedness of the path.
"""
CircularPolygon(p::AbstractPath;kw...) = CircularPolygon(ClosedPath(p;kw...))
function CircularPolygon(p::AbstractVector{T};kw...) where T<:AbstractCurve
	CircularPolygon(ClosedPath(p;kw...))
end

# Required methods
curves(p::CircularPolygon) = curves(p.path)
curve(p::CircularPolygon,k::Integer) = curve(p.path,k)
arclength(p::CircularPolygon) = arclength(p.path)
(p::CircularPolygon)(t) = point(p.path,t)

inv(p::CircularPolygon) = CircularPolygon([inv(c) for c in curves(p)])

# Other methods
"""
	ispositive(p::CircularPolygon)
Determine whether the circular polygon is positively oriented (i.e., circulates counterclockwise around the points it encloses).
"""
ispositive(p::CircularPolygon) = winding(0,1/p) < 0

#
# Polygon
#

# Type
"""
	(type) Polygon
Type for closed paths consisting entirely of segments and rays.
"""
struct Polygon <: AbstractPolygon
	path
	function Polygon(p::AbstractClosedPath)
		# Assumes continuity and closure have been checked previously
		valid = isa.(curves(p),Union{Segment,Ray})
		@assert all(valid) "All sides must be a Segment or Ray"
		new(p)
	end
end

# Constructors
"""
	Polygon(p::AbstractPath; tol=<default>)
	Polygon(p::AbstractVector{T<:AbstractCurve}; tol=<default>)
Construct a polygon from a (possibly closed) path, or from a vector of curves. The `tol` parameter is a tolerance used when checking continuity and closedness of the path.
"""
Polygon(p::AbstractPath;kw...) = Polygon(ClosedPath(p;kw...))
function Polygon(p::AbstractVector{T};kw...) where T<:AbstractCurve
	Polygon(ClosedPath(p;kw...))
end

"""
	Polygon(v::AbstractVector)
Construct a polygon from a vector of its vertices. Each element of `v` should be either a finite vertex, or a tuple of two angles that indicate the angles of two rays incident to an infinite vertex: one "to" infinity, and a second "from" infinity.
"""
function Polygon(v::AbstractVector)
	n = length(v)
	p = Vector{Union{Segment,Ray}}(undef,n)
	for j = 1:n
		vthis = v[j]
		vnext = v[mod(j,n)+1]
		if isa(vthis,Tuple)
			if isa(vnext,Tuple)
				@error("Cannot have consecutive infinite vertices")
				return nothing
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
	return Polygon(ClosedPath(p))
end

# Required methods
curves(p::Polygon) = curves(p.path)
arclength(p::Polygon) = arclength(p.path)
(p::Polygon)(t) = point(p.path,t)

inv(p::Polygon) = CircularPolygon([inv(c) for c in curves(p)])

# Display methods
function show(io::IO,::MIME"text/plain",P::Polygon)
	print(io,"Polygon with ",length(P)," vertices:")
	for (v,a) in zip(vertices(P),angles(P))
		print(io,"\n   ")
		show(io,MIME("text/plain"),v)
		print(io,", interior angle ",a/pi,"π")
	end
end

# Other methods
"""
	truncate(P::Union{CircularPolygon,Polygon})
Apply `truncate` to `P` using a circle that is centered at the centroid of its finite vertices, and a radius twice the maximum from the centroid to the finite vertices.
"""
function truncate(p::Union{CircularPolygon,Polygon})
	isfinite(p) && return p   # nothing to do
	return truncate(p,enclosing_circle(p,4))
end

"""
	truncate(P::Union{CircularPolygon,Polygon},C::Circle)
Compute a trucated form of the polygon by replacing each pair of rays incident at infinity with two segments connected by an arc along the given circle. This is *not* a true clipping of the polygon, as finite sides are not altered. The result is either a CircularPolygon or the original `P`.
"""
function truncate(p::Union{CircularPolygon,Polygon},c::Circle)
	n = length(p)
	s,v = sides(p),vertices(p)
	snew = Vector{Any}(undef,n)
	z_pre = NaN
	if isa(s[1],Ray) && isinf(point(s[1],0))
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
				t0 = arg(c,z_pre)
				delta = mod( angle( (z_post-c.center)/(z_pre-c.center) ) / (2π), 1)
				snew[k] = [Arc(c,t0,delta),Segment(z_post,s[k].base)]
				z_pre = NaN
			end
		end
	end
	return CircularPolygon(vcat(snew...))
end

"""
	angles(P::Polygon)
Compute a vector of interior angles at the vertices of the polygon `P`. At a finite vertex these lie in (0,2π]; at an infinite vertex, the angle is in [-2π,0].
"""
function angles(p::Polygon)
	# computes a turn angle in (-pi,pi]  (neg = left turn)
	turn(s1,s2) = π - mod(angle(s2/s1)+π,2π)
	s = unittangent.(p)
	n = length(p)
	v = vertices(p)
	θ = similar(real(s))
	for k = 1:n
		kprev = mod(k-2,n)+1
		θ[k] = π + turn(s[kprev],s[k])
		if isinf(v[k])
			θ[k] -= 2π
			if θ[k]==0
				# need a finite perturbation to distinguish 0,-2
				R = maximum(abs.(filter(isfinite,v)))
				C = Circle(0,100*R)
				zprev = intersect(p[kprev],C)
				znext = intersect(p[k],C)
				if angle(znext[1]/zprev[1]) < 0
					θ[k] -= 2π
				end
			end
		end
	end
	# # correct for possible clockwise orientation
	# sum(θ) > 0 ? θ : -θ
	return θ
end

"""
	ispositive(p::Polygon)
Determine whether the polygon is positively oriented (i.e., circulates counterclockwise around the points it encloses).
"""
ispositive(p::Polygon) = sum(angles(p)/pi .- 1) < 0

## Special polygon types and constructors

struct Rectangle <: AbstractPolygon
	center::Complex
	radii::SVector{2}
	rotation::Real   # angle about the center point
	polygon::Polygon
	function Rectangle(center::Number=0, radii::AbstractVector=[1.0, 1.0], rotation::Real=0)
		# rotation is in radians
		@assert (length(radii) == 2) && all(radii .>= 0)
		z = center .+ complex.(radii[1]*[-1,1,1,-1], radii[2]*[-1,-1,1,1])
		p = Polygon(z .* cis(rotation))
		new(complex(center), radii, rotation, p)
	end
end

function Base.show(io::IO,::MIME"text/plain", R::Rectangle)
	print(io,"Rectangle with center $(R.center), dimensions $(2R.radii[1]) × $(2R.radii[2]), rotation $(R.rotation/π)π")
end

function Base.show(io::IO, R::Rectangle)
	print(io,"$(2R.radii[1]) × $(2R.radii[2]) Rectangle")
end

# converters
Base.convert(::Type{Polygon}, r::Rectangle) = r.polygon
Polygon(r::Rectangle) = r.polygon

# required methods
arclength(r::Rectangle) = arclength(Polygon(r))
curves(r::Rectangle) = curves(Polygon(r))
point(r::Rectangle, t::Real) = point(Polygon(r), t)
(r::Rectangle)(t) = point(Polygon(r), t)
inv(r::Rectangle) = inv(Polygon(r))

# these provide shortcuts
Base.:+(r::Rectangle, z::Number) = Rectangle(r.center + z, r.radii, r.rotation)
Base.:+(z::Number, r::Rectangle) = r + z
Base.:-(r::Rectangle) = rectangle(-verices(r))
Base.:-(r::Rectangle, z::Number) = Rectangle(r.center - z, r.radii, r.rotation)
Base.:-(z::Number, r::Rectangle) = (-r) + z

# other methods
vertices(r::Rectangle) = vertices(Polygon(r))
angles(::Rectangle) = fill(π/2, 4)
function Base.extrema(r::Rectangle)
	z = vertices(r)
	return extrema(real(z)), extrema(imag(z))
end

# alternate constructors
"""
	rectangle(v)
Construct the rectangle with vertices given in the vector `v`.
"""
function rectangle(v::AbstractVector{<:Number})
	@assert length(v) == 4
	p = Polygon(v)
	if !ispositive(p)
		p = reverse(p)
		v = reverse(v)
	end
	@assert all( isapprox.(angles(p), π/2, atol=1e-13) ) "Given vertices do not form a rectangle"
	center = sum(v) / 4
	v .-= center
	θ = angle(v[2] - v[1])
	v *= cis(-θ)
	radii = abs.([real(v[2] - v[1]), imag(v[4] - v[1])])
	return Rectangle(center, radii, θ)
end

"""
	rectangle(xlim, ylim)
Construct the rectangle defined by `xlim[1]` < Re(z) < `xlim[2]`, `ylim[1]` < Im(z) < `ylim[2]`.
"""
function rectangle(xlim::AbstractVector, ylim::AbstractVector)
	@assert (length(xlim)==2) && (length(ylim)==2)
	@assert (xlim[1] < xlim[2]) && (ylim[1] < ylim[2])
	center = complex(mean(xlim), mean(ylim))
	radii = SVector((xlim[2] - xlim[1])/2, (ylim[2] - ylim[1])/2)
	return Rectangle(center, radii, 0)
	# x = [xlim[1],xlim[2],xlim[2],xlim[1],xlim[1]]
	# y = [ylim[1],ylim[1],ylim[2],ylim[2],ylim[1]]
	# Polygon( [Segment(complex(x[k],y[k]),complex(x[k+1],y[k+1])) for k in 1:4] )
end

"""
	rectangle(z1, z2)
Construct the axes-aligned rectangle whose opposing corners are the given complex values.
"""
rectangle(z1::AnyComplex, z2::AnyComplex) = rectangle(sort(real([z1,z2])), sort(imag([z1,z2])))
rectangle(z1::Number, z2::Number) = rectangle(promote(complex(float(z1)), complex(float(z2)))...)

"""
	n_gon(n)
Construct a regular n-gon with vertices on the unit circle.
"""
function n_gon(n::Integer)
	@assert n > 2 "Must have at least three vertices"
	Polygon( exp.(2im*pi*(0:n-1)/n) )
end
