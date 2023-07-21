const CR = ComplexRegions
const GB = Makie.GeometryBasics
export complex_theme
import .Makie
using ColorSchemes
using .Makie: PointBased, Poly, Lines, Series, Point2f, Combined, with_theme
import .Makie: convert_arguments, plottype, plot!

complex_theme = Makie.Theme(
    Axis = (aspect = Makie.DataAspect(),),
    Series = (linewidth = 4, color = ColorSchemes.seaborn_colorblind[1:10]),
    Lines = (color = ColorSchemes.seaborn_colorblind[1], linewidth=4),
    patchcolor = ColorSchemes.seaborn_colorblind[1],  # for regions
    Poly = (strokecolor=:black, strokewidth=4),       # also for regions
    )

# Allow plot of any complex vector
z_to_point(z::AnyComplex) = Point2f(reim(z)...)
convert_arguments(::PointBased, z::AbstractVector{<:Complex}) = (z_to_point.(z), )

#####
##### Curves and paths
#####

# Convert a pathlike object to a vector of points
plottype(::AbstractCurveOrPath) = Lines
plottype(::AbstractVector{<:AbstractCurveOrPath}) = Series
curve_to_points(c::AbstractCurveOrPath) = z_to_point.(plotdata(c))
convert_arguments(::PointBased, c::AbstractCurveOrPath) = (curve_to_points(c), )

# Plot a compound boundary (i.e., from a generic ConnectedRegion)
Compound = Tuple{Union{Nothing,AbstractJordan}, Vector{<:AbstractJordan}}
# plottype(::Compound) = Series
function plot!(plt::Combined{Any, S} where S<:Tuple{Compound})
    outer, inner = plt[1][]
    if !isnothing(outer)
        Makie.plot!(plt, outer)
    end
    Makie.plot!(plt, inner)
    return plt
end

#####
##### Regions
#####

# Convert a generic region to a Makie Polygon
plottype(::AbstractRegion) = Poly
function convert_arguments(PT::Type{<:Poly}, R::ConnectedRegion{N}) where N
    @show PT
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
