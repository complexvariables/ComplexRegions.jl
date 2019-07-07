module ComplexRegions

export Polar,Spherical
using ComplexValues,LinearAlgebra

AllComplex(T::Type{S}) where {S<:AbstractFloat} = Union{Complex{T},Polar{T},Spherical{T}}
AnyComplex{S} = Union{Complex{S},Polar{S},Spherical{S}} where {S<:AbstractFloat}

import Base: +,-,*,/,sign,inv,angle,real,imag,conj,show,iterate,eltype,length,getindex,isapprox,intersect

export inf
include("utilities.jl")

export point,start,stop,arclength,slope,isapprox,dist,closest,isleft,plotdata
export Curve,ClosedCurve,Circle,Line,Arc,Segment,Ray
include("curves.jl")

export Path,ClosedPath,Polygon,CircularPolygon
export curve,vertex,side,breakindex,isbounded,sign,angle,winding
include("paths.jl")

include("plotrecipes.jl")

end # module
