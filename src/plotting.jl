const CR = ComplexRegions
const GB = Makie.GeometryBasics
using ColorSchemes, Infiltrator
using Makie: PointBased, convert_arguments, plottype

complex_theme = Makie.Theme(
    Axis = (aspect = Makie.DataAspect(),),
    linewidth = 6,
    colormap = ColorSchemes.seaborn_colorblind,
    palette = (patchcolor = ColorSchemes.seaborn_colorblind[1:10],),
    patchcolor = ColorSchemes.seaborn_colorblind[1],
    Poly = (strokecolor=:black, strokewidth=6),
    )

# Allow plot of any complex vector
z_to_point(z::AnyComplex) = Makie.Point2f(reim(z)...)
Makie.convert_arguments(::PointBased, z::AbstractVector{<:Complex}) = (z_to_point.(z), )

# Convert a pathlike object to a vector of points
plottype(::AbstractCurveOrPath) = Makie.Lines
curve_to_points(c::AbstractCurveOrPath) = z_to_point.(plotdata(c))
convert_arguments(::PointBased, c::AbstractCurveOrPath) = (curve_to_points(c), )

# Convert a generic region to a Makie Polygon
plottype(::AbstractRegion) = Makie.Poly
function convert_arguments(PT::Type{<:Makie.Poly}, R::ConnectedRegion{N}) where N
    outer = curve_to_points(outerboundary(R))
    inner = curve_to_points.(innerboundary(R))
    return convert_arguments(PT, GB.Polygon(outer, inner))
end

# Conversions for particular cases
function convert_arguments(PT::Type{<:Makie.Poly}, R::InteriorSimplyConnectedRegion)
    ∂R = curve_to_points(boundary(R))
    return convert_arguments(PT, GB.Polygon(∂R))
end

function convert_arguments(PT::Type{<:Makie.Poly}, R::ExteriorSimplyConnectedRegion)
    return convert_arguments(PT, truncate(R))
end

function convert_arguments(PT::Type{<:Makie.Poly}, R::ExteriorRegion{N}) where N
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

convert_arguments(PT::Type{<:Makie.Poly}, A::Annulus) = convert_arguments(PT, ConnectedRegion{2}(A.outer, [A.inner]))

# RegionPlot = Makie.Combined{Any, S} where {N, S <: Tuple{ConnectedRegion{N}}}
# function Makie.plot!(plt::RegionPlot)
#     pol = Makie.convert_single_argument(plt[1][])
#     Makie.poly!(plt, pol)
#     plt
# end

# ISCPlot = Makie.Combined{Any, S} where {S <: Tuple{InteriorSimplyConnectedRegion}}
# function Makie.plot!(plt::ISCPlot)
#     pol = Makie.convert_single_argument(plt[1][])
#     Makie.poly!(plt, pol)
#     plt
# end

# ESCPlot = Makie.Combined{Any, S} where {S <: Tuple{ExteriorSimplyConnectedRegion}}
# function Makie.plot!(plt::ESCPlot)
#     R = truncate(plt[1][])
#     Makie.plot!(plt, R)
#     C = outerboundary(R)
#     zc = C.center
#     r = 0.15*C.radius
#     if !isnothing(Makie.current_axis())
#         Makie.current_axis().limits = (real(zc) - r, real(zc) + r, imag(zc) - r, imag(zc) + r)
#     end
#     plt
# end

# AnnPlot = Makie.Combined{Any, S} where {N, S <: Tuple{Annulus}}
# function Makie.plot!(plt::AnnPlot)
#     pol = Makie.convert_single_argument(plt[1][])
#     Makie.plot!(plt, pol)
#     plt
# end
