AbstractJordan = Union{AbstractClosedCurve,AbstractClosedPath}
abstract type AbstractRegion end

""" 
	isfinite(R::AbstractRegion) 

Return `true` if the region is bounded in the complex plane.
"""
isfinite(r::AbstractRegion) = @error "No isfinite() method defined for type $(typeof(r))"

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

#
# ConnectedRegion 
#

abstract type AbstractConnectedRegion{N} <: AbstractRegion end
AbstractSimplyConnectedRegion = AbstractConnectedRegion{1}

struct SimplyConnectedRegion{T<:AbstractJordan} <: AbstractConnectedRegion{1}
	boundary::T 
end

boundary(R::SimplyConnectedRegion) = R.boundary
isfinite(R::SimplyConnectedRegion) = isfinite(boundary(R))

function show(io::IO,R::SimplyConnectedRegion)
	print(IOContext(io,:compact=>true),"Region to the left of ",R.boundary)
end
function show(io::IO,::MIME"text/plain",R::SimplyConnectedRegion)
	print(io,"Region to the left of:\n   ",R.boundary)
end

in(z::Number,R::SimplyConnectedRegion) = isleft(z,R.boundary)
!(R::SimplyConnectedRegion) = SimplyConnectedRegion(reverse(R.boundary))
isapprox(R1::SimplyConnectedRegion,R2::SimplyConnectedRegion) = R1.boundary ≈ R2.boundary

struct ConnectedRegion{N} <: AbstractConnectedRegion{N}
	outer::Union{Nothing,AbstractJordan} 
	inner::AbstractVector 
	function ConnectedRegion{N}(outer,inner) where N
		if isnothing(outer)
			@assert N == length(inner) 
		else 
			@assert N == length(inner) + 1
		end 
		new(outer,inner)
	end
end

function in(z::Number,R::ConnectedRegion) 
	val = all( !isleft(z,c) for c in R.inner )
	isnothing(R.outer) ? val : (val && isleft(z,R.outer))
end

function boundary(R::ConnectedRegion)
	isnothing(R.outer) ? R.inner : R.outer,R.inner 
end

#
# special cases
#

region(C::AbstractJordan,left=true) = SimplyConnectedRegion{typeof(C)}(C)
interior(C::AbstractJordan) = region(C)
exterior(C::AbstractJordan) = region(reverse(C))
between(outer::AbstractJordan,inner::AbstractJordan) = ConnectedRegion{2}(outer,inner)

AbstractDisk = SimplyConnectedRegion{T} where T<:Circle
disk(C::Circle) = interior(C) 
disk(center::Number,radius::Real) = interior(Circle(center,radius))
unitdisk = disk(complex(0.0),1.0)
function show(io::IO,::MIME"text/plain",R::AbstractDisk)
	side = in(Inf,R) ? "exterior" : "interior"
	print(io,"Disk $side to:\n   ",R.boundary)
end

AbstractHalfplane = SimplyConnectedRegion{T} where T<:Line 
halfplane(L::Line) = interior(L)
halfplane(a::Number,b::Number) = interior(Line(a,b))
upperhalfplane = halfplane(Line(0.0,direction=1.0))
lowerhalfplane = halfplane(Line(0.0,direction=-1.0))
lefthalfplane = halfplane(Line(0.0,direction=1.0im))
righthalfplane = halfplane(Line(0.0,direction=-1.0im))
function show(io::IO,::MIME"text/plain",R::AbstractHalfplane)
	print(io,"Half-plane to the left of:\n   ",R.boundary)
end

PolygonalRegion = SimplyConnectedRegion{Polygon} 

struct Annulus{S,T} <: AbstractConnectedRegion{2} 
	outer::Circle{S} 
	inner::Circle{T} 
end
function Annulus(center::Number,outerrad::Real,innerrad::Real)  
	@assert outerrad > innerrad > 0 
	Annulus(Circle(center,outerrad),Circle(center,innerrad))
end
Annulus(outerrad::Real,innerrad::Real) = Annulus(complex(0.0),outerrad,innerrad)

in(z::Number,A::Annulus) = isleft(z,A.outer) && !isleft(z,A.inner)
boundary(A::Annulus) = A.outer,A.inner 
isfinite(A::Annulus) = isfinite(A.outer) 

function show(io::IO,::MIME"text/plain",R::Annulus)
	print(io,"Annulus interior to:\n   ",R.outer,"\nand exterior to:\n   ",R.inner)
end
