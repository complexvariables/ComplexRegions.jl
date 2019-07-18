module ComplexRegions

export Polar,Spherical
using ComplexValues,LinearAlgebra,StaticArrays

AllComplex(T::Type{S}) where {S<:AbstractFloat} = Union{Complex{T},Polar{T},Spherical{T}}
AnyComplex{S} = Union{Complex{S},Polar{S},Spherical{S}} where {S<:AbstractFloat}

import Base: +,-,*,/,!,∘,sign,inv,angle,real,imag,conj,show,iterate,eltype,length,getindex,isapprox,intersect,union,truncate,reverse,in

export inf
include("utilities.jl")

export point,arclength,slope,isapprox,dist,closest,isleft,plotdata,isright,reflect,tangent,normal,arg
export Curve,ClosedCurve,Circle,Line,Arc,Segment,Ray
include("curves.jl")

export Path,ClosedPath,Polygon,CircularPolygon
export curve,vertex,side,breakindex,isbounded,sign,angle,winding,rectangle
include("paths.jl")

export SimplyConnectedRegion,ConnectedRegion,disk,unitdisk,Annulus,PolygonalRegion 
export halfplane,upperhalfplane,lowerhalfplane,lefthalfplane,righthalfplane
export RegionIntersection,RegionUnion
export region,interior,exterior,between,boundary
include("regions.jl")

export Möbius,Mobius
include("mobius.jl")

include("plotrecipes.jl")

end # module
