#
# generic curves
#

abstract type AbstractCurve end

point(c::AbstractCurve,t::Real) = @error "No point() method defined for type $(typeof(c))"
point(c::AbstractCurve,t::AbstractArray{T}) where T<:Real = [point(c,t) for t in t]
start(c::AbstractCurve) = point(c,0.0)
stop(c::AbstractCurve) = point(c,1.0)
arclength(c::AbstractCurve) = @error "No arclength() method defined for type $(typeof(c))"
plotdata(C::AbstractCurve) = point(C,LinRange(0,1,501))

struct Curve <: AbstractCurve 
	point 
	arclength 
end
Curve(f) = Curve(f,missing)
Curve(f,a::Real,b::Real,arclen=missing) = Curve(t -> f(scaleto(a,b,t)),arclen)

point(C::Curve,t::Real) = C.point(t)
(C::Curve)(t::Real) = point(C,t)
arclength(C::Curve) = C.arclength 
conj(C::Curve) = Curve(t->conj(C.point(t)),C.arclength)
reverse(C::Curve) = Curve(t->C.point(1-t),C.arclength)

#
# closed curves
# 

abstract type AbstractClosedCurve <: AbstractCurve end

struct ClosedCurve <: AbstractClosedCurve 
	point 
	arclength 
	function ClosedCurve(f,arclen=missing;tol=1e-12)
		@assert isapprox(f(0),f(1);rtol=tol,atol=tol) "Curve does not close"
		new(f,arclen)
	end
end
ClosedCurve(f,a::Real,b::Real,arclen=missing;kw...) = ClosedCurve(t -> f(scaleto(a,b,t)),arclen;kw...)

point(C::ClosedCurve,t::Real) = C.point(t)
(C::ClosedCurve)(t::Real) = point(C,t)
arclength(C::ClosedCurve) = C.arclength 
conj(C::ClosedCurve) = ClosedCurve(t->conj(C.point(t)),C.arclength)
reverse(C::ClosedCurve) = ClosedCurve(t->C.point(1-t),C.arclength)

include("lines_segments.jl")
include("circles_arcs.jl")

# Intersections
function twolines_meet(z1,s1,z2,s2,tol) 
	M = [ real(s1) -real(s2); imag(s1) -imag(s2) ]
	if cond(M) > 1/tol
		# parallel lines, caller must decide 
		return NaN,NaN
	else
		d = z2 - z1
		t = M \ [real(d);imag(d)]
		return t[1],t[2]
	end
end

function intersect(l1::Line,l2::Line;tol=1e-12) 
	s1,s2 = l1.direction,l2.direction
	z1,z2 = l1.base,l2.base
	t1,t2 = twolines_meet(z1,s1,z2,s2,tol)
	if isnan(t1)     # parallel lines
		# either identical or no intersection
		return dist(z2,l1) < tol*(abs(z2)+1) ? l1 : []
	else
		return z1 + t1*s1 
	end
end

function intersect(g1::Segment,g2::Segment;tol=1e-12)
	z1,z2 = g1.za,g2.za 
	s1,s2 = g1.zb-z1,g2.zb-z2 
	t1,t2 = twolines_meet(z1,s1,z2,s2,tol)
	if isnan(t1)   # parallel lines
		# scale first segment to [0,1]
		a = real( (g2.za-g1.za)/s1 )
		b = real( (g2.zb-g1.za)/s1 )
		if a > b 
			c = b; b = a; a = c 
		end
		if b < -tol || a > 1+tol
			return []
		else
			return Segment(z1+max(0,a)*s1,z1+min(1,b)*s1)
		end
	else   # nonparallel
		if (-tol ≤ t1 ≤ 1+tol) && (-tol ≤ t2 ≤ 1+tol)
			return z1 + t1*s1 
		else
			return [] 
		end
	end
end

intersect(g::Segment,l::Line;tol=1e-12) = intersect(l,g,tol=tol)
function intersect(l::Line,g::Segment;tol=1e-12)
	z1,z2 = l.base,g.za
	s1,s2 = l.direction,g.zb-g.za
	t1,t2 = twolines_meet(z1,s1,z2,s2,tol)
	if isnan(t1)   # parallel lines
		return dist(g.za,l) < tol*(1+abs(g.za)) ? g : [] 
	else
		return 0 ≤ t2 ≤ 1 ? z1+t1*s1 : []
	end
end

function intersect(r1::Ray,r2::Ray;tol=1e-12) 
	s1,s2 = exp(complex(0,r1.angle)),exp(complex(0,r2.angle))
	z1,z2 = r1.base,r2.base
	t1,t2 = twolines_meet(z1,s1,z2,s2,tol)
	if isnan(t1)     # parallel lines
		if dist(z2,r1) < tol*(abs(z2)+1)
			return r2 
		else
			dist(z1,r2) < tol*(abs(z1)+1) ? r1 : []
		end
	else
		return (t1 ≥ 0) && (t2 ≥ 0) ? z1 + t1*s1 : []
	end
end

intersect(r::Ray,l::Line;tol=1e-12) = intersect(l,r,tol=tol)
function intersect(l::Line,r::Ray;tol=1e-12)
	z1,z2 = l.base,r.base
	s1,s2 = l.direction,exp(1im*r.angle)
	t1,t2 = twolines_meet(z1,s1,z2,s2,tol)
	if isnan(t1)   # parallel lines
		return dist(r.base,l) < tol*(1+abs(r.base)) ? r : [] 
	else
		return 0 ≤ t2 ? z1+t1*s1 : []
	end
end

intersect(g::Segment,r::Ray;tol=1e-12) = intersect(r,g,tol=tol)
function intersect(r::Ray,g::Segment;tol=1e-12)
	z1,z2 = r.base,g.za
	s1,s2 = exp(1im*r.base),g.zb-g.za
	t1,t2 = twolines_meet(z1,s1,z2,s2,tol)
	if isnan(t1)   # parallel lines
		# move ray to positive Re axis
		a = real((g.za-z1)/s1)
		b = real((g.zb-z1)/s1)
		if a > b 
			c = b; b = a; a = c 
		end
		return b < -tol ? [] : Segment(z1+max(0,a)*s1,z1+b*s1)
	else   # nonparallel
		return (-tol ≤ t1) && (-tol ≤ t2 ≤ 1+tol) ? z1 + t1*s1  : []
	end
end


function intersect(c1::Circle,c2::Circle;tol=1e-12) 
	r1,r2 = c1.radius,c2.radius 
	z1,z2 = c1.center,c2.center
	delta = z2-z1
	d = abs(delta) 
	if abs(d) < tol 
		if isapprox(r1,r2,rtol=tol,atol=tol) 
			return c1 
		else
			return []
		end
	elseif (d > r1+r2) || d < abs(r1-r2) 
		return [] 
	else
		a = (r1^2 - r2^2 + d^2)/(2d) 
		p = z1 + a*delta/d 
		h = sqrt(r1^2-a^2) 
		w = 1im*h*delta/d
		return [p+w,p-w]
	end
end

intersect(l::Line,c::Circle;tol=1e-12) = intersect(c,l,tol=tol)
function intersect(c::Circle,l::Line;tol=1e-12)
	# find a radius that is perpendicular to the line
	p = Line(c.center,direction=1im*l.direction)
	z = intersect(p,l) 
	a,r = abs(z-c.center),c.radius
	if a ≤ r + tol 
		c = sqrt(r^2-a^2)*l.direction
		return [z+c,z-c] 
	else 
		return []
	end
end

intersect(r::Ray,c::Circle;tol=1e-12) = intersect(c,r,tol=tol)
function intersect(c::Circle,r::Ray;tol=1e-12)
	# probably not the most efficient way
	z = intersect(c,Line(r.base,direction=exp(1im*r.angle)),tol=tol)
	d = [ dist(z,r) for z in z ]
	return z[d.<=tol]
end

intersect(s::Segment,c::Circle;tol=1e-12) = intersect(c,s,tol=tol)
function intersect(c::Circle,s::Segment;tol=1e-12)
	# probably not the most efficient way
	z = intersect(c,Line(s.za,s.zb),tol=tol)
	d = [ dist(z,s) for z in z ]
	return z[d.<=tol]
end
