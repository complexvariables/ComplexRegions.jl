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
