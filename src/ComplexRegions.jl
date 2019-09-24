module ComplexRegions

DEFAULT = Dict(:tol=>1e-12)
"""
	ComplexRegions.default()
Return a dictionary of global default settings for the ComplexRegions package.

	ComplexRegions.default(key=value)
Change a global default setting in the running instance of the ComplexRegions package.
"""
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

using Reexport
@reexport using ComplexValues
using LinearAlgebra,StaticArrays

AnyComplex = Union{Complex{S},Polar{S},Spherical{S}} where {S<:AbstractFloat}

export Polar,Spherical
import Base: +,-,*,/,!,∘,sign,inv,angle,real,imag,conj,show,iterate,eltype,length,getindex,isapprox,isfinite,intersect,union,truncate,reverse,in
export +,-,*,/,!,∘,sign,inv,angle,real,imag,conj,show,iterate,eltype,length,getindex,isapprox,isfinite,intersect,union,truncate,reverse,in

include("utilities.jl")

export point,arclength,slope,dist,closest,isleft,plotdata,isright,reflect,tangent,unittangent,normal,arg,isinside,isoutside
export Curve,ClosedCurve,Circle,Line,Arc,Segment,Ray
include("curves.jl")

export Path,ClosedPath,Polygon,CircularPolygon
export curve,curves,vertex,vertices,side,sides,isfinite,sign,angle,angles,winding,ispositive,rectangle,n_gon
include("paths.jl")

export SimplyConnectedRegion,ConnectedRegion,ExteriorRegion,disk,unitdisk,Annulus,PolygonalRegion 
export halfplane,upperhalfplane,lowerhalfplane,lefthalfplane,righthalfplane
export RegionIntersection,RegionUnion
export region,interior,exterior,between,boundary,innerboundary,outerboundary
include("regions.jl")
include("simplyconnected.jl")

export Möbius,Mobius
include("mobius.jl")

include("docs.jl")
include("plotrecipes.jl")

end # module
