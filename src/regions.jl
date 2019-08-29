AbstractJordan = Union{AbstractClosedCurve,AbstractClosedPath}
abstract type AbstractRegion end

# Required methods
"""
	boundary(R::AbstractRegion)
Return the boundary of a region. Depending on the type of region, this might be a vector.
"""
boundary(R::AbstractRegion) = @error "No boundary() method defined for type $(typeof(R))"

"""
	in(z::Number,R::AbstractRegion;tol=<default>)
	z ∈ R   (type "\\in" followed by tab)
True if `z` is in the region `R`. 
"""
in(z::Number,R::AbstractRegion;tol=DEFAULT[:tol]) = @error "No in() method defined for type $(typeof(R))"
	
# Default implementations
""" 
	isfinite(R::AbstractRegion) 
Return `true` if the region is bounded in the complex plane.
"""
isfinite(R::AbstractRegion) = @error "No isfinite() method defined for type $(typeof(R))"

"""
	(type) RegionIntersection 
Representation of the intersection of two regions.
"""
struct RegionIntersection <: AbstractRegion
	one::AbstractRegion
	two::AbstractRegion 
end
in(z::Number,R::RegionIntersection) = in(z,R.one) && in(z,R.two)

"""
	(type) RegionUnion 
Representation of the union of two regions.
"""
struct RegionUnion <: AbstractRegion
	one::AbstractRegion
	two::AbstractRegion 
end
in(z::Number,R::RegionUnion) = in(z,R.one) || in(z,R.two)

"""
	intersect(R1::AbstractRegion,R2::AbstractRegion)
	R1 ∩ R2    (type "\\cap" followed by tab key)
Create the region that is the intersection of `R1` and `R2`. 
"""
intersect(R1::AbstractRegion,R2::AbstractRegion) = RegionIntersection(R1,R2)

"""
	union(R1::AbstractRegion,R2::AbstractRegion)
	R1 ∪ R2    (type "\\cup" followed by tab key)
Create the region that is the union of `R1` and `R2`. 
"""
union(R1::AbstractRegion,R2::AbstractRegion) = RegionUnion(R1,R2)

# 
# AbstractConnectedRegion
#

abstract type AbstractConnectedRegion{N} <: AbstractRegion end

# Required methods
innerboundary(R::AbstractConnectedRegion) = @error "No innerboundary() method defined for type $(typeof(R))"
outerboundary(R::AbstractConnectedRegion) = @error "No outerboundary() method defined for type $(typeof(R))"
boundary(R::AbstractConnectedRegion) = outerboundary(R),innerboundary(R)

#
# SimplyConnectedRegion 
#

AbstractSimplyConnectedRegion = AbstractConnectedRegion{1}

"""
	(type) SimplyConnectedRegion 
Representation of a simply connected region in the extended complex plane. 
	SimplyConnectedRegion(p::Union{AbstractClosedCurve,AbstractClosedPath})
Construct an open simply connected region by specifying its boundary. The region is "to the left" of the orientation of the boundary.
"""
struct InteriorSimplyConnectedRegion{T<:AbstractJordan} <: AbstractConnectedRegion{1}
	boundary::T 
end

struct ExteriorSimplyConnectedRegion{T<:AbstractJordan} <: AbstractConnectedRegion{1}
	boundary::T 
end

SimplyConnectedRegion = Union{InteriorSimplyConnectedRegion{T},ExteriorSimplyConnectedRegion{T}} where T<:AbstractJordan

"""
	interior(C)
Construct the region interior to the closed curve or path `C`. If `C` is bounded, the bounded enclosure is chosen regardless of the orientation of `C`; otherwise, the region "to the left" is the interior. 
"""
function interior(C::AbstractJordan) 
	if isfinite(C) && winding(1/C,0) > 0
		C = reverse(C)
	end
	InteriorSimplyConnectedRegion(C)
end

"""
	exterior(C)
Construct the region exterior to  the closed curve or path `C`. If `C` is bounded, the bounded enclosure is chosen regardless of the orientation of `C`; otherwise, the region "to the right" is the exterior. 
"""
function exterior(C::AbstractJordan) 
	if isfinite(C) 
		if winding(1/C,0) < 0
			C = reverse(C)
		end
	else
		if C isa AbstractClosedPath && (length(filter(isinf,vertices(C))) > 1)
			@error "Disconnected exterior"
		end
		C = reverse(C)
	end
	ExteriorSimplyConnectedRegion(C)
end

boundary(R::SimplyConnectedRegion) = R.boundary
innerboundary(R::InteriorSimplyConnectedRegion) = nothing
innerboundary(R::ExteriorSimplyConnectedRegion) = R.boundary
outerboundary(R::InteriorSimplyConnectedRegion) = R.boundary
outerboundary(R::ExteriorSimplyConnectedRegion) = nothing
in(z::Number,R::InteriorSimplyConnectedRegion) = isinside(z,outerboundary(R))
in(z::Number,R::ExteriorSimplyConnectedRegion) = isoutside(z,innerboundary(R))
isfinite(R::InteriorSimplyConnectedRegion) = isfinite(outerboundary(R))
isfinite(R::ExteriorSimplyConnectedRegion) = false

function show(io::IO,R::InteriorSimplyConnectedRegion)
	print(IOContext(io,:compact=>true),"Region interior to ",R.boundary)
end

function show(io::IO,R::ExteriorSimplyConnectedRegion)
	print(IOContext(io,:compact=>true),"Region exterior to ",R.boundary)
end

function show(io::IO,::MIME"text/plain",R::SimplyConnectedRegion)
	show(io,R)
end

"""
	!(R::SimplyConnectedRegion)
Compute the region complementary to `R`. This is not quite set complementation, as neither region includes its boundary. The complement is always simply connected in the extended plane. 
"""
!(R::InteriorSimplyConnectedRegion) = exterior(reverse(boundary(R)))
!(R::ExteriorSimplyConnectedRegion) = interior(reverse(boundary(R)))

"""
	isapprox(R1::SimplyConnectedRegion,R2::SimplyConnectedRegion; tol=<default>)
Determine whether `R1` and `R2` represent the same region, up to tolerance `tol`. Equivalently, determine whether their boundaries are the same.
"""
isapprox(R1::S,R2::T;tol=DEFAULT[:tol]) where {S,T<:SimplyConnectedRegion} = isapprox(boundary(R1),boundary(R2),tol=tol)

#
# ExteriorRegion
#
struct ExteriorRegion{N} <: AbstractConnectedRegion{N}
	inner::AbstractVector 
	function ExteriorRegion{N}(inner) where N
		@assert N == length(inner) "Incorrect connectivity"
		@assert all(c isa AbstractJordan for c in inner) "Boundary components must be closed curves or paths"
		@assert all(isfinite.(inner)) "Inner boundaries must be finite"
		# correct orientations of inner components
		b = copy(inner)
		for c in b 
			if isoutside(Inf,c)
				 c = reverse(c)
			end
		end
		new(b)
	end
end  

ExteriorRegion(inner) = ExteriorRegion{length(inner)}(inner)
in(z::Number,R::ExteriorRegion) = all( isoutside(z,c) for c in R.inner )
innerboundary(R::ExteriorRegion) = R.inner 
outerboundary(R::ExteriorRegion) = nothing 

#
# ConnectedRegion 
#

"""
	(type) ConnectedRegion{N} 
Representation of a `N`-connected region in the extended complex plane. 
"""
struct ConnectedRegion{N} <: AbstractConnectedRegion{N}
	outer::Union{Nothing,AbstractJordan} 
	inner::AbstractVector 
	function ConnectedRegion{N}(outer,inner) where N
		n = length(inner) + !isnothing(outer)
		@assert N == n "Incorrect connectivity"
		if !isnothing(outer)
			# correct orientation of outer component 
			isin = [isinside(point(c,0),outer) for c in inner ]
			if all(.!isin)
				outer = reverse(outer) 
			else
				@assert all(isin) "Inner components appear to be crossing the outer boundary"
			end 
		end
		# correct orientations of inner components
		@assert all(isfinite.(inner)) "Inner boundaries must be finite"
		for c in inner 
			if !isfinite(interior(c))
				c = reverse(c) 
			end
		end
		new(outer,inner)
	end
end  

"""
	ConnectedRegion(outer,inner)
Construct an open connected region by specifying its boundary components. The `outer` boundary could be `nothing` or a closed curve or path. The `inner` boundary should be a vector of one or more nonintersecting closed curves or paths. The defined region is interior to the outer boundary and exterior to all the components of the inner boundary, regardless of the orientations of the given curves. 
"""
function ConnectedRegion(outer,inner) 
	n = length(inner)
	if isnothing(outer) || isempty(outer) 
		ExteriorRegion{n}(inner)
	else
		ConnectedRegion{n+1}(outer,inner)
	end
end

function in(z::Number,R::ConnectedRegion) 
	all( isoutside(z,c) for c in R.inner ) && isinside(z,R.outer)
end

outerboundary(R::ConnectedRegion) = R.outer
innerboundary(R::ConnectedRegion) = R.inner 
innerboundary(R::ConnectedRegion{2}) = R.inner[1] 

#
# special cases
#

"""
	between(outer,inner)
Construct the region interior to the closed curve or path `outer` and interior to `inner`. 
"""
function between(outer::AbstractJordan,inner::AbstractJordan)
	if isfinite(outer) && isinside(Inf,outer)
		outer = reverse(outer)
	end
	if isfinite(inner) && isoutside(Inf,inner)
		inner = reverse(inner)
	end
	ConnectedRegion{2}(outer,[inner])
end

# disks
AbstractDisk = SimplyConnectedRegion{T} where T<:Circle
"""
	disk(C::Circle) 
Construct the disk interior to `C`.
"""
disk(C::Circle) = interior(C) 
"""
	disk(center::Number,radius::Real) 
Construct the disk with the given `center` and `radius`. 
"""
disk(center::Number,radius::Real) = interior(Circle(center,radius))
unitdisk = disk(complex(0.0),1.0)
function show(io::IO,::MIME"text/plain",R::AbstractDisk)
	side = in(Inf,R) ? "exterior" : "interior"
	print(io,"Disk $side to:\n   ",R.boundary)
end

# half-planes
AbstractHalfplane = SimplyConnectedRegion{T} where T<:Line
"""
	halfplane(L::Line) 
Construct the half-plane to the left of `L`.
"""
halfplane(L::Line) = interior(L)
"""
	halfplane(a,b) 
Construct the half-plane to the left of the line from `a` to `b`.
"""
halfplane(a::Number,b::Number) = interior(Line(a,b))
upperhalfplane = halfplane(Line(0.0,direction=1.0))
lowerhalfplane = halfplane(Line(0.0,direction=-1.0))
lefthalfplane = halfplane(Line(0.0,direction=1.0im))
righthalfplane = halfplane(Line(0.0,direction=-1.0im))
function show(io::IO,::MIME"text/plain",R::AbstractHalfplane)
	print(io,"Half-plane to the left of:\n   ",R.boundary)
end

"""
	(type) PolygonalRegion
Representation of a simply connected region bounded by a plane.
"""
PolygonalRegion = SimplyConnectedRegion{Polygon} 

# Annulus
"""
	(type) Annulus 
Representation of the region between two circles.
"""
struct Annulus <: AbstractConnectedRegion{2} 
	outer::Circle
	inner::Circle
	function Annulus(outer::Circle,inner::Circle) 
		@assert(outer.center ≈ inner.center)
		if isinside(Inf,outer)
			outer = reverse(outer)
		end
		if isoutside(Inf,inner)
			inner = reverse(inner)
		end
		new(outer,inner)
	end
end
"""
	Annulus(radouter,radinner)
	Annulus(radouter,radinner,center)
Construct a concentric annulus of outer radius `radouter` and inner radius `radinner` centered at `center`. If the center is not given, the origin is used.
"""
function Annulus(outerrad::Real,innerrad::Real,center::Number=0)  
	@assert outerrad > innerrad > 0 
	Annulus(Circle(center,outerrad,true),Circle(center,innerrad,false))
end

innerboundary(A::Annulus) = A.inner 
outerboundary(A::Annulus) = A.outer
in(z::Number,A::Annulus) = isinside(z,A.outer) && isoutside(z,A.inner)
isfinite(::Annulus) = true

function show(io::IO,::MIME"text/plain",R::Annulus)
	print(io,"Annulus in the complex plane:\n")
	print(io,"   centered at ",R.outer.center," with distances from ",R.inner.radius," to ",R.outer.radius)
end
