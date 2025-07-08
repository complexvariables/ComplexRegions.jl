const AbstractSimplyConnectedRegion{T,S} = AbstractConnectedRegion{1,T}

struct InteriorSimplyConnectedRegion{T,S} <: AbstractSimplyConnectedRegion{T} where S<:AbstractJordan{T}
    boundary::S
end

struct ExteriorSimplyConnectedRegion{T,S} <: AbstractSimplyConnectedRegion{T} where S<:AbstractJordan{T}
    boundary::S
end

"""
	(type) SimplyConnectedRegion
Representation of a simply connected region in the extended complex plane.
"""
const SimplyConnectedRegion{T,S} = Union{InteriorSimplyConnectedRegion{T,S},ExteriorSimplyConnectedRegion{T,S}}

"""
	interior(C)
Construct the region interior to the closed curve or path `C`. If `C` is bounded, the bounded enclosure is chosen regardless of the orientation of `C`; otherwise, the region "to the left" is the interior.
"""
function interior(C::AbstractJordan{T}) where {T}
    # Determine the winding number of the inverse path around the origin. This reveals the
    # orientation without needing to find an interior point.
    # use a strange number to avoid divide by zero when it's on the curve
    if isfinite(C) && winding(1 / (C + complex(T(0.13298), -T(0.398127))), 0) > 0
        C = reverse(C)
    end
    InteriorSimplyConnectedRegion{T,typeof(C)}(C)
end

"""
	exterior(C)
Construct the region exterior to  the closed curve or path `C`. If `C` is bounded, the bounded enclosure is chosen regardless of the orientation of `C`; otherwise, the region "to the right" is the exterior.
"""
function exterior(C::AbstractJordan{T}) where {T}
    if isfinite(C)
        # use a strange number to avoid divide by zero when it's on the curve
        if winding(1 / (C + complex(T(0.13298), -T(0.398127))), 0) < 0
            C = reverse(C)
        end
    else
        if C isa AbstractClosedPath && (length(filter(isinf, vertices(C))) > 1)
            @error "Disconnected exterior"
        end
        #C = reverse(C)
    end
    ExteriorSimplyConnectedRegion{T,typeof(C)}(C)
end

boundary(R::SimplyConnectedRegion) = R.boundary
innerboundary(R::InteriorSimplyConnectedRegion) = nothing
innerboundary(R::ExteriorSimplyConnectedRegion) = R.boundary
outerboundary(R::InteriorSimplyConnectedRegion) = R.boundary
outerboundary(R::ExteriorSimplyConnectedRegion) = nothing
in(z::Number, R::InteriorSimplyConnectedRegion) = isinside(z, outerboundary(R))
in(z::Number, R::ExteriorSimplyConnectedRegion) = isoutside(z, innerboundary(R))
isfinite(R::InteriorSimplyConnectedRegion) = isfinite(outerboundary(R))
isfinite(R::ExteriorSimplyConnectedRegion) = false

function show(io::IO, R::InteriorSimplyConnectedRegion)
    print(IOContext(io, :compact => true), "Region interior to ", R.boundary)
end

function show(io::IO, R::ExteriorSimplyConnectedRegion)
    print(IOContext(io, :compact => true), "Region exterior to ", R.boundary)
end

function show(io::IO, ::MIME"text/plain", R::SimplyConnectedRegion)
    show(io, R)
end

"""
	!(R::SimplyConnectedRegion)
Compute the region complementary to `R`. This is not quite set complementation, as neither region includes its boundary. The complement is always simply connected in the extended plane.
"""
function Base.:!(R::InteriorSimplyConnectedRegion{T,S}) where {T,S}
	return ExteriorSimplyConnectedRegion{T,S}(boundary(R))
end

function Base.:!(R::ExteriorSimplyConnectedRegion{T,S}) where {T,S}
	return InteriorSimplyConnectedRegion{T,S}(boundary(R))
end

"""
	isapprox(R1::SimplyConnectedRegion, R2::SimplyConnectedRegion; tol=<default>)
Determine whether `R1` and `R2` represent the same region, up to tolerance `tol`. Equivalently, determine whether their boundaries are the same.
"""
function isapprox(
				R1::SimplyConnectedRegion{T,U},
				R2::SimplyConnectedRegion{S,R};
				tol=tolerance(S,T)
				) where {S,T,U,R}
	return isapprox(boundary(R1), boundary(R2); tol=tol)
end

# disks
const AbstractDisk{T} = InteriorSimplyConnectedRegion{T,Circle{T}}
"""
	disk(C::Circle)
Construct the disk interior to `C`.
"""
disk(C::Circle) = interior(C)
"""
	disk(center::Number,radius::Real)
Construct the disk with the given `center` and `radius`.
"""
disk(center::Number, radius::Real) = interior(Circle(center, radius))
unitdisk = disk(complex(0.0), 1.0)

# COV_EXCL_START
function show(io::IO, ::MIME"text/plain", R::AbstractDisk)
    side = in(Inf, R) ? "exterior" : "interior"
    print(io, "Disk $side to:\n   ", R.boundary)
end
# COV_EXCL_END

# quads
AbstractQuad{T} = InteriorSimplyConnectedRegion{T,Rectangle{T}}
"""
	quad(R::Rectangle)
Construct the rectangle interior to `R`.
"""
quad(R::Rectangle) = interior(R)
unitquad = quad(Rectangle(0.0, [1, 1]))

# COV_EXCL_START
function show(io::IO, ::MIME"text/plain", R::AbstractQuad)
    print(io, "Quad inside ", R.boundary)
end
# COV_EXCL_END

# half-planes
AbstractHalfplane{T} = SimplyConnectedRegion{T,Line{T}}
"""
	halfplane(L::Line)
Construct the half-plane to the left of `L`.
"""
halfplane(L::Line) = interior(L)
"""
	halfplane(a,b)
Construct the half-plane to the left of the line from `a` to `b`.
"""
halfplane(a::Number, b::Number) = interior(Line(a, b))
upperhalfplane = halfplane(Line(0.0, direction=1.0))
lowerhalfplane = halfplane(Line(0.0, direction=-1.0))
lefthalfplane = halfplane(Line(0.0, direction=1.0im))
righthalfplane = halfplane(Line(0.0, direction=-1.0im))
function show(io::IO, ::MIME"text/plain", R::AbstractHalfplane)
    print(io, "Half-plane to the left of:\n   ", R.boundary)
end

"""
	(type) PolygonalRegion
Representation of a simply connected region bounded by a polygon.
"""
PolygonalRegion{T} = SimplyConnectedRegion{T,Polygon{T}}
