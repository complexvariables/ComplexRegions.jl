AbstractJordan = Union{AbstractClosedCurve,AbstractClosedPath}
abstract type AbstractRegion end

# Required methods
boundary(R::AbstractRegion) = @error "No boundary() method defined for type $(typeof(R))"

"""
	in(z::Number,R::AbstractRegion;tol=<default>)
	z ∈ R   (type "\\in" followed by tab)
True if `z` is in the region `R`. 
"""
in(z::Number,R::AbstractRegion;tol=DEFAULT[:tol]) = @error "No in() method defined for type $(typeof(R))"
	
""" 
	isfinite(R::AbstractRegion) 
Return `true` if the region is bounded in the complex plane.
"""
isfinite(R::AbstractRegion) = @error "No isfinite() method defined for type $(typeof(R))"

# Default implementations

"""
	(type) RegionIntersection 
Representation of the intersection of two regions.
"""
struct RegionIntersection <: AbstractRegion
	one::AbstractRegion
	two::AbstractRegion 
end
in(z::Number,R::RegionIntersection) = in(z,R.one) && in(z,R.two)

"""
	(type) RegionUnion 
Representation of the union of two regions.
"""
struct RegionUnion <: AbstractRegion
	one::AbstractRegion
	two::AbstractRegion 
end
in(z::Number,R::RegionUnion) = in(z,R.one) || in(z,R.two)

"""
	intersect(R1::AbstractRegion,R2::AbstractRegion)
	R1 ∩ R2    (type "\\cap" followed by tab key)
Create the region that is the intersection of `R1` and `R2`. 
"""
intersect(R1::AbstractRegion,R2::AbstractRegion) = RegionIntersection(R1,R2)

"""
	union(R1::AbstractRegion,R2::AbstractRegion)
	R1 ∪ R2    (type "\\cup" followed by tab key)
Create the region that is the union of `R1` and `R2`. 
"""
union(R1::AbstractRegion,R2::AbstractRegion) = RegionUnion(R1,R2)

# 
# AbstractConnectedRegion
#

abstract type AbstractConnectedRegion{N} <: AbstractRegion end

# Required methods
innerboundary(R::AbstractConnectedRegion) = @error "No innerboundary() method defined for type $(typeof(R))"
outerboundary(R::AbstractConnectedRegion) = @error "No outerboundary() method defined for type $(typeof(R))"
boundary(R::AbstractConnectedRegion) = outerboundary(R),innerboundary(R)

# Default implementations

+(R::AbstractConnectedRegion,z::Number) = typeof(R)(outerboundary(R)+z,innerboundary(R).+z)
+(z::Number,R::AbstractConnectedRegion) = +(R,z)

-(R::AbstractConnectedRegion) = typeof(R)(-outerboundary(R),-innerboundary(R))
-(R::AbstractConnectedRegion,z::Number) = +(R,-z)
-(z::Number,R::AbstractConnectedRegion) = +(z,-R)

*(R::AbstractConnectedRegion,z::Number) = typeof(R)(outerboundary(R)*z,innerboundary(R)*z)
*(z::Number,R::AbstractConnectedRegion) = R*z

/(R::AbstractConnectedRegion,z::Number) = *(R,1/z)
#/(z::Number,R::AbstractConnectedRegion) = z*inv(R)
#inv(p::AbstractConnectedRegion) = typeof(R)([inv(c) for c in curves(p)])

#
# concrete implementations
#

#
# ExteriorRegion
#
struct ExteriorRegion{N} <: AbstractConnectedRegion{N}
	inner::AbstractVector 
	function ExteriorRegion{N}(inner) where N
		@assert N == length(inner) "Incorrect connectivity"
		@assert all(c isa AbstractJordan for c in inner) "Boundary components must be closed curves or paths"
		@assert all(isfinite.(inner)) "Inner boundaries must be finite"
		# correct orientations of inner components
		b = copy(inner)
		for c in b 
			if isoutside(Inf,c)
				 c = reverse(c)
			end
		end
		new(b)
	end
end  

ExteriorRegion(inner) = ExteriorRegion{length(inner)}(inner)
in(z::Number,R::ExteriorRegion) = all( isoutside(z,c) for c in R.inner )
innerboundary(R::ExteriorRegion) = R.inner 
outerboundary(R::ExteriorRegion) = nothing 

#
# ConnectedRegion 
#

"""
	(type) ConnectedRegion{N} 
Representation of a `N`-connected region in the extended complex plane. 
"""
struct ConnectedRegion{N} <: AbstractConnectedRegion{N}
	outer::Union{Nothing,AbstractJordan} 
	inner::AbstractVector 
	function ConnectedRegion{N}(outer,inner) where N
		n = length(inner) + !isnothing(outer)
		@assert N == n "Incorrect connectivity"
		if !isnothing(outer)
			# correct orientation of outer component?
			isin = [isinside(point(c,0),outer) for c in inner ]
			if all(.!isin)
				outer = reverse(outer) 
			else
				@assert all(isin) "Inner components appear to be crossing the outer boundary"
			end 
		end
		# correct orientations of inner components?
		@assert all(isfinite.(inner)) "Inner boundaries must be finite"
		for c in inner 
			if !isfinite(interior(c))
				c = reverse(c) 
			end
		end
		new(outer,inner)
	end
end  

"""
	ConnectedRegion(outer,inner)
Construct an open connected region by specifying its boundary components. The `outer` boundary could be `nothing` or a closed curve or path. The `inner` boundary should be a vector of one or more nonintersecting closed curves or paths. The defined region is interior to the outer boundary and exterior to all the components of the inner boundary, regardless of the orientations of the given curves. 
"""
function ConnectedRegion(outer,inner) 
	n = length(inner)
	if isnothing(outer) || isempty(outer) 
		ExteriorRegion{n}(inner)
	else
		ConnectedRegion{n+1}(outer,inner)
	end
end

function in(z::Number,R::ConnectedRegion) 
	all( isoutside(z,c) for c in R.inner ) && isinside(z,R.outer)
end

outerboundary(R::ConnectedRegion) = R.outer
innerboundary(R::ConnectedRegion) = R.inner 
#innerboundary(R::ConnectedRegion{2}) = R.inner[1] 

#
# special cases
#

"""
	between(outer,inner)
Construct the region interior to the closed curve or path `outer` and interior to `inner`. 
"""
function between(outer::AbstractJordan,inner::AbstractJordan)
	if isfinite(outer) && isinside(Inf,outer)
		outer = reverse(outer)
	end
	if isfinite(inner) && isoutside(Inf,inner)
		inner = reverse(inner)
	end
	ConnectedRegion{2}(outer,[inner])
end


# Annulus

"""
	(type) Annulus 
Representation of the region between two circles.
"""
struct Annulus <: AbstractConnectedRegion{2} 
	outer::Circle
	inner::Circle
	function Annulus(outer::Circle,inner::Circle) 
		@assert(outer.center ≈ inner.center)
		if isinside(Inf,outer)
			outer = reverse(outer)
		end
		if isoutside(Inf,inner)
			inner = reverse(inner)
		end
		new(outer,inner)
	end
end
"""
	Annulus(radouter,radinner)
	Annulus(radouter,radinner,center)
Construct a concentric annulus of outer radius `radouter` and inner radius `radinner` centered at `center`. If the center is not given, the origin is used.
"""
function Annulus(outerrad::Real,innerrad::Real,center::Number=0)  
	@assert outerrad > innerrad > 0 
	Annulus(Circle(center,outerrad,true),Circle(center,innerrad,false))
end

innerboundary(A::Annulus) = A.inner 
outerboundary(A::Annulus) = A.outer
in(z::Number,A::Annulus) = isinside(z,A.outer) && isoutside(z,A.inner)
isfinite(::Annulus) = true

function show(io::IO,::MIME"text/plain",R::Annulus)
	print(io,"Annulus in the complex plane:\n")
	print(io,"   centered at ",R.outer.center," with distances from ",R.inner.radius," to ",R.outer.radius)
end
