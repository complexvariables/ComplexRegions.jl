@doc """
	intersect(c1::AbstractCurve,c2::AbstractCurve; tol=<default>)
Find the intersection(s) of two curves. The result could be a vector of zero or more values, or a curve. 
""" intersect

# utility, to be used in what follows
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

function intersect(l1::Line,l2::Line;tol=DEFAULT[:tol]) 
	s1,s2 = l1.direction,l2.direction
	z1,z2 = l1.base,l2.base
	t1,t2 = twolines_meet(z1,s1,z2,s2,tol)
	if isnan(t1)     # parallel lines
		# either identical or no intersection
		return dist(z2,l1) < tol*(abs(z2)+1) ? l1 : []
	else
		return [z1 + t1*s1]
	end
end

function intersect(g1::Segment,g2::Segment;tol=DEFAULT[:tol])
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
			return [z1 + t1*s1]
		else
			return [] 
		end
	end
end

intersect(g::Segment,l::Line;tol=DEFAULT[:tol]) = intersect(l,g,tol=tol)
function intersect(l::Line,g::Segment;tol=DEFAULT[:tol])
	z1,z2 = l.base,g.za
	s1,s2 = l.direction,g.zb-g.za
	t1,t2 = twolines_meet(z1,s1,z2,s2,tol)
	if isnan(t1)   # parallel lines
		return dist(g.za,l) < tol*(1+abs(g.za)) ? g : []
	else
		return 0 ≤ t2 ≤ 1 ? [z1+t1*s1] : []
	end
end

function intersect(r1::Ray,r2::Ray;tol=DEFAULT[:tol]) 
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
		return (t1 ≥ 0) && (t2 ≥ 0) ? [z1 + t1*s1] : []
	end
end

intersect(r::Ray,l::Line;tol=DEFAULT[:tol]) = intersect(l,r,tol=tol)
function intersect(l::Line,r::Ray;tol=DEFAULT[:tol])
	z1,z2 = l.base,r.base
	s1,s2 = l.direction,exp(1im*r.angle)
	t1,t2 = twolines_meet(z1,s1,z2,s2,tol)
	if isnan(t1)   # parallel lines
		return dist(r.base,l) < tol*(1+abs(r.base)) ? r : [] 
	else
		return 0 ≤ t2 ? [z1+t1*s1] : []
	end
end

intersect(g::Segment,r::Ray;kw...) = intersect(r,g;kw...)
function intersect(r::Ray,g::Segment;tol=DEFAULT[:tol])
	z1,z2 = r.base,g.za
	s1,s2 = exp(1im*r.angle),g.zb-g.za
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
		return (-tol ≤ t1) && (-tol ≤ t2 ≤ 1+tol) ? [z1 + t1*s1]  : []
	end
end


function intersect(c1::Circle,c2::Circle;tol=DEFAULT[:tol]) 
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
		z = [p+w,p-w]
		return z
	end
end

intersect(l::Line,c::Circle;kw...) = intersect(c,l;kw...)
function intersect(c::Circle,l::Line;tol=DEFAULT[:tol])
	# find a radius that is perpendicular to the line
	p = Line(c.center,direction=1im*l.direction)
	zi = intersect(p,l)[1] 
	a,r = abs(zi-c.center),c.radius
	if a ≤ r + tol 
		c = sqrt(r^2-a^2)*l.direction
		if abs(c) > tol
			return [zi+c,zi-c]  
		else
			return [zi]
		end 
	else 
		return []
	end
end

intersect(r::Ray,c::Circle;tol=DEFAULT[:tol]) = intersect(c,r,tol=tol)
function intersect(c::Circle,r::Ray;tol=DEFAULT[:tol])
	# probably not the most efficient way
	z = intersect(c,Line(r.base,direction=exp(1im*r.angle)),tol=tol)
	d = [ dist(z,r) for z in z ]
	return z[d.<=tol]
end

intersect(s::Segment,c::Circle;tol=DEFAULT[:tol]) = intersect(c,s,tol=tol)
function intersect(c::Circle,s::Segment;tol=DEFAULT[:tol])
	# probably not the most efficient way
	z = intersect(c,Line(s.za,s.zb),tol=tol)
	d = [ dist(z,s) for z in z ]
	return z[d.<=tol]
end

function intersect(a1::Arc,a2::Arc;tol=DEFAULT[:tol]) 
	if isapprox(a1.circle,a2.circle,tol=tol)
		# to what extent do the arcs overlap on the same circle? 
		r = a1.circle.radius
		d10ison2 = dist(point(a1,0),a2) < tol*(1+r)
		d11ison2 = dist(point(a1,1),a2) < tol*(1+r)
		d20ison1 = dist(point(a2,0),a1) < tol*(1+r)
		d21ison1 = dist(point(a2,1),a1) < tol*(1+r)
		if d10ison2 
			if d11ison2
				return a1 
			else
				t = d20ison1 ? arg(a1,point(a2,0)) : arg(a1,point(a2,1))
				return Arc(point(a1,[0,t/2,t]))
			end 
		elseif d11ison2 
			t = d20ison1 ? arg(a1,point(a2,0)) : arg(a1,point(a2,1))
			return Arc(point(a1,[1-t,1-t/2,1]))
		else
			return [] 
		end 
	else 	
		z = intersect(a1.circle,a2.circle) 
		f = w -> (dist(w,a1) < tol) & (dist(w,a2) < tol)
		return filter(f,z)
	end
end 

intersect(c::AbstractCurve,a::Arc;kw...) = intersect(a,c;kw...)
function intersect(a::Arc,c::AbstractCurve;tol=DEFAULT[:tol]) 
	if isapprox(a.circle,c,tol=tol)
		return a 
	else
		z = intersect(a.circle,c)
		return filter(v->dist(v,a)<tol,z)
	end
end

# Determine the (directional) crossings of a Ray(y*1im,0) with a given curve.
# For intersection at an endpoint, give half the value, which gives the right winding number for a polygon.
function raycrossing(z,c::Segment;tol=DEFAULT[:tol])
	s = 0
	y = imag(z) 
	a,b = imag(c(0)),imag(c(1))
	atend = abs(a-y) < tol || abs(b-y) < tol
	if y ≥ a
		if y ≤ b && isleft(z,c)
			s = atend ? s+0.5 : s+1
		end
	else
		if y ≥ b && !isleft(z,c) 
			s = atend ? s-0.5 : s-1
		end
	end
	return s
end

function raycrossing(z,c::AbstractCurve;tol=DEFAULT[:tol])
	s = 0
	for z in intersect(c,Ray(z,0)) 
		t = arg(c,z) 
		if -tol ≤ t ≤ 1+tol
			y = imag(unittangent(c,t))
			if abs(y) < tol 
				# horiztonal tangent; check a perturbed location
				ϵ = sqrt(tol)
				if t < 0.5 
					y = imag(unittangent(c,t+2ϵ))
				else 
					y = -imag(unittangent(c,t-2ϵ))
				end
			end
			d = tol < t < 1-tol ? copysign(1,y) : copysign(1,y)/2
			s += d
		end
	end
	return s
end
