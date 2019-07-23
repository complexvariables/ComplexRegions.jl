module ComplexRegions

DEFAULT = Dict(:tol=>1e-12)
function default(;kw...) 
	if isempty(kw) 
		return DEFAULT
	end
	for (sym,val) in kw
		if !haskey(DEFAULT,sym)
			@error "Unrecognized default setting `$sym`" 
		else
			DEFAULT[sym] = val
			@info "Default value of `$sym` set to $val."
		end
	end
	return nothing
end

export Polar,Spherical
using ComplexValues,LinearAlgebra,StaticArrays

AllComplex(T::Type{S}) where {S<:AbstractFloat} = Union{Complex{T},Polar{T},Spherical{T}}
AnyComplex{S} = Union{Complex{S},Polar{S},Spherical{S}} where {S<:AbstractFloat}

import Base: +,-,*,/,!,∘,sign,inv,angle,real,imag,conj,show,iterate,eltype,length,getindex,isapprox,isfinite,intersect,union,truncate,reverse,in

include("utilities.jl")

export point,arclength,slope,isapprox,dist,closest,isleft,plotdata,isright,reflect,tangent,normal,arg
export Curve,ClosedCurve,Circle,Line,Arc,Segment,Ray
include("curves.jl")

export Path,ClosedPath,Polygon,CircularPolygon
export curve,vertex,side,isfinite,sign,angle,winding,rectangle
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
