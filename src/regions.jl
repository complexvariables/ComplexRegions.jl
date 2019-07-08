AbstractJordan = Union{AbstractClosedCurve,AbstractClosedPath}
abstract type AbstractRegion end

struct Region{T<:AbstractJordan} <: AbstractRegion
	boundary::T 
	left::Bool 
end
Region(C::AbstractJordan,left=true) = Region{typeof(C)}(C,left)

in(z::Number,R::Region) = !xor(R.left,isleft(z,R.boundary))

Disk = Region{Circle{T}} where T<:AnyComplex
Halfplane = Region{Line{T}} where T<:AnyComplex
PolygonalRegion = Region{Polygon} 

struct RegionIntersection <: AbstractRegion
	one::AbstractRegion
	two::AbstractRegion 
end
in(z::Number,R::RegionIntersection) = in(z,R.one) && in(z,R.two)

struct RegionUnion <: AbstractRegion
	one::AbstractRegion
	two::AbstractRegion 
end
in(z::Number,R::RegionUnion) = in(z,R.one) || in(z,R.two)

intersect(R1::AbstractRegion,R2::AbstractRegion) = RegionIntersection(R1,R2)
union(R1::AbstractRegion,R2::AbstractRegion) = RegionUnion(R1,R2)
