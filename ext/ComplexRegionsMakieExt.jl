module ComplexRegionsMakieExt
using Makie
using ComplexRegions
const GB = Makie.GeometryBasics

const AbstractCurve = ComplexRegions.AbstractCurve
const AbstractPath = ComplexRegions.AbstractPath
const AbstractCurveOrPath = Union{AbstractCurve, AbstractPath}
const AbstractJordan = ComplexRegions.AbstractJordan
const AbstractCircularPolygon = ComplexRegions.AbstractCircularPolygon
const AbstractRegion = ComplexRegions.AbstractRegion
const ExteriorRegion = ComplexRegions.ExteriorRegion
const ExteriorSimplyConnectedRegion = ComplexRegions.ExteriorSimplyConnectedRegion
const InteriorSimplyConnectedRegion = ComplexRegions.InteriorSimplyConnectedRegion

z_to_point(z::Complex{T} where T) = Makie.Point2f(reim(z)...)

# Allow plot of any complex vector
Makie.convert_arguments(::PointBased, z::AbstractVector{<:Complex}) = (z_to_point.(z), )

##### Curves and paths

# Convert a pathlike object to a vector of points
Makie.plottype(::AbstractCurveOrPath) = Lines
Makie.plottype(::AbstractVector{<:AbstractCurveOrPath}) = Series
curve_to_points(c::AbstractCurveOrPath) = z_to_point.(complex(plotdata(c)))
Makie.convert_arguments(::PointBased, c::AbstractCurveOrPath) = (curve_to_points(c), )

# Plot a compound boundary (e.g., from a generic ConnectedRegion)
const Compound = Tuple{Union{Nothing,AbstractJordan}, Vector{<:AbstractJordan}}
function Makie.plot!(plt::Tuple{Compound})
    outer, inner = plt[1][]
    if !isnothing(outer)
        Makie.plot!(plt, outer)
    end
    Makie.plot!(plt, inner)
    return plt
end

##### Regions

# Convert a generic region to a Makie Polygon
Makie.plottype(::AbstractRegion) = Poly
function Makie.convert_arguments(PT::Type{<:Poly}, R::InteriorConnectedRegion{N}) where N
    outer = curve_to_points(outerboundary(R))
    inner = curve_to_points.(innerboundary(R))
    return convert_arguments(PT, GB.Polygon(outer, inner))
end

# Conversions for particular cases
function Makie.convert_arguments(PT::Type{<:Poly}, R::InteriorSimplyConnectedRegion)
    ∂R = curve_to_points(boundary(R))
    return convert_arguments(PT, GB.Polygon(∂R))
end

function Makie.convert_arguments(PT::Type{<:Poly}, R::ExteriorSimplyConnectedRegion)
    return convert_arguments(PT, truncate(R))
end

function Makie.convert_arguments(PT::Type{<:Poly}, R::ExteriorRegion{N}) where N
    return convert_arguments(PT, truncate(R))
end

function Base.truncate(R::ExteriorSimplyConnectedRegion)
    ∂R = boundary(R)
    C = ComplexRegions.enclosing_circle(ClosedPath(∂R), 8)
    return between(reverse(C), ∂R)
end

function Base.truncate(R::ExteriorRegion)
    ∂R = innerboundary(R)
    C = ComplexRegions.enclosing_circle(ClosedPath.(∂R), 8)
    return connected_region(C, ∂R)
end

Makie.convert_arguments(PT::Type{<:Poly}, A::Annulus) = convert_arguments(PT, InteriorConnectedRegion{2}(A.outer, [A.inner]))

end
