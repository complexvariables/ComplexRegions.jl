module ComplexRegions

using Statistics
using Dierckx
using CircularArrays

using Reexport
@reexport using ComplexValues
using LinearAlgebra, StaticArrays, ForwardDiff

AnyComplex{S} = Union{Complex{S},Polar{S},Spherical{S}} where {S<:AbstractFloat}

tolerance(S::Type{<:AbstractFloat}, T::Type{<:AbstractFloat}=S) = 100max(eps(S),eps(T))

import ComplexValues: real_type
# The following is more convenient than convert(Polar{T}, ::Polar{S}), etc.,
# because it doesn't require specifying the complex type
convert_real_type(T::Type{<:Real}, z::Complex{S}) where S = Complex{T}(z)
convert_real_type(T::Type{<:Real}, z::Polar{S}) where S = Polar{T}(z)
convert_real_type(T::Type{<:Real}, z::Spherical{S}) where S = Spherical{T}(z)
convert_real_type(T::Type{<:Real}, x::S) where S<:Number = convert(T, x)

import Base: !, ∘, sign, inv, angle, real, imag, conj, show, iterate, eltype, length, getindex, isapprox, isfinite, intersect, union, truncate, reverse, in
export !, ∘, sign, inv, angle, real, imag, conj, show, iterate, eltype, length, getindex, isapprox, isfinite, intersect, union, truncate, reverse, in

include("utilities.jl")

export point, arclength, slope, dist, closest, isleft, isright, reflect, tangent, unittangent, normal, arg, isinside, isoutside, isclosed, plotdata, unitcircle
export Curve, ClosedCurve, Circle, Line, Arc, Segment, Ray
include("curves.jl")

export Path, ClosedPath, Polygon, polygon, CircularPolygon
export curve, curves, vertex, vertices, side, sides, isfinite, sign, angle, angles, winding, ispositive, Rectangle, rectangle, n_gon
include("paths.jl")

const AbstractCurveOrPath{T} = Union{AbstractCurve{T},AbstractPath{T}}

export SimplyConnectedRegion, ConnectedRegion, ExteriorRegion, disk, unitdisk, Annulus, PolygonalRegion
export halfplane, upperhalfplane, lowerhalfplane, lefthalfplane, righthalfplane
export RegionIntersection, RegionUnion
export region, interior, exterior, between, boundary, innerboundary, outerboundary, modulus
include("regions.jl")
include("simplyconnected.jl")

export discretize
include("discretize.jl")

export Möbius,Mobius
include("mobius.jl")

export Shapes
include("shapes.jl")

include("docs.jl")

end # module
