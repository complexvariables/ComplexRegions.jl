AbstractJordan = Union{AbstractClosedCurve,AbstractClosedPath}
abstract type AbstractRegion end
abstract type AbstractConnectedRegion{N} <: AbstractRegion end

""" 
	isfinite(R::AbstractRegion) 
Return `true` if the region is bounded in the complex plane.
"""
isfinite(r::AbstractRegion) = @error "No isfinite() method defined for type $(typeof(r))"

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
# SimplyConnectedRegion 
#

AbstractSimplyConnectedRegion = AbstractConnectedRegion{1}

"""
	(type) SimplyConnectedRegion 
Representation of a simply connected region in the extended complex plane. 
	SimplyConnectedRegion(p::Union{AbstractClosedCurve,AbstractClosedPath})
Construct an open simply connected region by specifying its boundary. The region is "to the left" of the orientation of the boundary.
"""
struct SimplyConnectedRegion{T<:AbstractJordan} <: AbstractConnectedRegion{1}
	boundary::T 
end

"""
	boundary(R::SimplyConnectedRegion)
Return the boundary of a simply connected region.
"""
boundary(R::SimplyConnectedRegion) = R.boundary
isfinite(R::SimplyConnectedRegion) = isfinite(boundary(R)) && winding(Inf,boundary(R))==0

function show(io::IO,R::SimplyConnectedRegion)
	print(IOContext(io,:compact=>true),"Region to the left of ",R.boundary)
end
function show(io::IO,::MIME"text/plain",R::SimplyConnectedRegion)
	print(io,"Region to the left of:\n   ",R.boundary)
end

"""
	in(z,R::SimplyConnectedRegion)
	z ∈ R    (type "\\in" followed by tab key)
Determine whether the given value `z` is in the region `R`.
"""
in(z::Number,R::SimplyConnectedRegion) = isleft(z,R.boundary)

"""
	!(R::SimplyConnectedRegion)
Compute the region complementary to `R`. This is not quite set complementation, as neither region includes its boundary. The complement is always simply connected in the extended plane. 
"""
!(R::SimplyConnectedRegion) = SimplyConnectedRegion(reverse(R.boundary))

"""
	isapprox(R1::SimplyConnectedRegion,R2::SimplyConnectedRegion; tol=<default>)
Determine whether `R1` and `R2` represent the same region, up to tolerance `tol`. Equivalently, determine whether their boundaries are the same.
"""
isapprox(R1::SimplyConnectedRegion,R2::SimplyConnectedRegion;tol=DEFAULT[:tol]) = isapprox(R1.boundary,R2.boundary,tol=tol)

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
			R = interior(outer)
			isin = [point(c(0)) ∈ R for c in inner ]
			if all(.!isin)
				outer = reverse(outer) 
			else
				@assert !all(isin) "Inner components appear to be crossing the outer boundary"
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
	n = length(inner) + !isnothing(outer)
	ConnectedRegion{n}(outer,inner)
end

function in(z::Number,R::ConnectedRegion) 
	val = all( !isleft(z,c) for c in R.inner )
	isnothing(R.outer) ? val : (val && isleft(z,R.outer))
end

boundary(R::ConnectedRegion) = R.outer,R.inner 

#
# special cases
#

region(C::AbstractJordan,left=true) = SimplyConnectedRegion{typeof(C)}(C)
"""
	interior(C)
Construct the region interior to (to the left of) the closed curve or path `C`. 
"""
interior(C::AbstractJordan) = region(C)
"""
	exterior(C)
Construct the region exterior to (to the right of) the closed curve or path `C`. 
"""
exterior(C::AbstractJordan) = region(reverse(C))
"""
	between(outer,inner)
Construct the region interior to the closed curve or path `outer` and interior to `inner`. 
"""
between(outer::AbstractJordan,inner::AbstractJordan) = ConnectedRegion{2}(outer,inner)

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
struct Annulus{S,T} <: AbstractConnectedRegion{2} 
	outer::Circle{S} 
	inner::Circle{T} 
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

in(z::Number,A::Annulus) = isleft(z,A.outer) && !isleft(z,A.inner)
boundary(A::Annulus) = A.outer,A.inner 
isfinite(A::Annulus) = isfinite(A.outer) 

function show(io::IO,::MIME"text/plain",R::Annulus)
	print(io,"Annulus interior to:\n   ",R.outer,"\nand exterior to:\n   ",R.inner)
end
