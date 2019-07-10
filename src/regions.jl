AbstractJordan = Union{AbstractClosedCurve,AbstractClosedPath}
abstract type AbstractRegion end

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
	left::Bool
	SimplyConnectedRegion{T}(C::T,left=true) where T<:AbstractJordan = new(C,left)
end
SimplyConnectedRegion(C::AbstractJordan,left=true) = SimplyConnectedRegion{typeof(C)}(C,left)

boundary(R::SimplyConnectedRegion) = R.boundary

function show(io::IO,R::SimplyConnectedRegion)
	dir = R.left ? "left" : "right"
	print(IOContext(io,:compact=>true),"Region to the $dir of ",R.boundary)
end
function show(io::IO,::MIME"text/plain",R::SimplyConnectedRegion)
	dir = R.left ? "left" : "right"
	print(io,"Region to the $dir of:\n   ",R.boundary)
end

in(z::Number,R::SimplyConnectedRegion) = !xor(R.left,isleft(z,R.boundary))

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

region(C::AbstractJordan,left=true) = SimplyConnectedRegion{typeof(C)}(C,left)
interior(C::AbstractJordan) = region(C,true)
exterior(C::AbstractJordan) = region(C,false)
between(outer::AbstractJordan,inner::AbstractJordan) = ConnectedRegion{2}(outer,inner)

Disk = SimplyConnectedRegion{Circle} 
SimplyConnectedRegion{Circle}(center::Number,radius::Real) = Disk(Circle(center,radius))
UnitDisk = Disk(complex(0.0),1.0)
function show(io::IO,::MIME"text/plain",R::Disk)
	side = R.left ? "interior" : "exterior"
	print(io,"Disk $side to:\n   ",R.boundary)
end

Halfplane = SimplyConnectedRegion{Line} 
UpperHalfplane = Halfplane(Line(0.0,direction=1.0))
LowerHalfplane = Halfplane(Line(0.0,direction=-1.0))
LeftHalfplane = Halfplane(Line(0.0,direction=1.0im))
RightHalfplane = Halfplane(Line(0.0,direction=-1.0im))
function show(io::IO,::MIME"text/plain",R::Halfplane)
	dir = R.left ? "left" : "right"
	print(io,"Half-plane to the $dir of:\n   ",R.boundary)
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

function show(io::IO,::MIME"text/plain",R::Annulus)
	print(io,"Annulus interior to:\n   ",R.outer,"\nand exterior to:\n   ",R.inner)
end
