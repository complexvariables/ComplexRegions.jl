const CR = ComplexRegions
const GB = Makie.GeometryBasics
export complex_theme
import .Makie
using ColorSchemes
using .Makie: PointBased, Poly, Lines, Point2f
import .Makie: convert_arguments, plottype

complex_theme = Makie.Theme(
    Axis = (aspect = Makie.DataAspect(),),
    linewidth = 5,
    colormap = ColorSchemes.seaborn_colorblind,
    palette = (patchcolor = ColorSchemes.seaborn_colorblind[1:10],),
    patchcolor = ColorSchemes.seaborn_colorblind[1],
    Poly = (strokecolor=:black, strokewidth=6),
    )

# Allow plot of any complex vector
z_to_point(z::AnyComplex) = Point2f(reim(z)...)
convert_arguments(::PointBased, z::AbstractVector{<:Complex}) = (z_to_point.(z), )

# Convert a pathlike object to a vector of points
plottype(::AbstractCurveOrPath) = Lines
curve_to_points(c::AbstractCurveOrPath) = z_to_point.(plotdata(c))
convert_arguments(::PointBased, c::AbstractCurveOrPath) = (curve_to_points(c), )

# Convert a generic region to a Makie Polygon
plottype(::AbstractRegion) = Poly
function convert_arguments(PT::Type{<:Poly}, R::ConnectedRegion{N}) where N
    outer = curve_to_points(outerboundary(R))
    inner = curve_to_points.(innerboundary(R))
    return convert_arguments(PT, GB.Polygon(outer, inner))
end

# Conversions for particular cases
function convert_arguments(PT::Type{<:Poly}, R::InteriorSimplyConnectedRegion)
    ∂R = curve_to_points(boundary(R))
    return convert_arguments(PT, GB.Polygon(∂R))
end

function convert_arguments(PT::Type{<:Poly}, R::ExteriorSimplyConnectedRegion)
    return convert_arguments(PT, truncate(R))
end

function convert_arguments(PT::Type{<:Poly}, R::ExteriorRegion{N}) where N
    return convert_arguments(PT, truncate(R))
end

function Base.truncate(R::ExteriorSimplyConnectedRegion)
    ∂R = boundary(R)
    C = enclosing_circle(ClosedPath(∂R), 8)
    return between(reverse(C), ∂R)
end

function Base.truncate(R::ExteriorRegion)
    ∂R = innerboundary(R)
    C = enclosing_circle(ClosedPath.(∂R), 8)
    return ConnectedRegion(C, ∂R)
end

convert_arguments(PT::Type{<:Poly}, A::Annulus) = convert_arguments(PT, ConnectedRegion{2}(A.outer, [A.inner]))
