module ComplexRegions

export Polar,Spherical
using ComplexValues,LinearAlgebra

AllComplex(T::Type{S}) where {S<:AbstractFloat} = Union{Complex{T},Polar{T},Spherical{T}}
AnyComplex = Union{Complex,Polar,Spherical}

import Base: +,-,*,/,sign,inv,angle,real,imag,conj,show

include("utilities.jl")

export point,start,stop,arclength,Circle,Segment
include("curves.jl")

end # module
