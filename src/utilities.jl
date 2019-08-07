# Scale from/to [0,1] 
scalefrom(a,b,t) = @. (t-a)/(b-a)
scaleto(a,b,t) = @. a + t*(b-a)

# Unique real roots of a quadratic.
function realroots(a,b,c)
	a==0 && return [-c/b]
	bh = -b/2
	d = bh^2 - a*c 
	if d < 0
		[]
	elseif d==0
		[bh/a]
	else
		r = (bh + sign(bh)*sqrt(d)) / a 
		[r,c/(r*a)]
	end
end  

# Are three points arranged in counterclockwise order? 
function isccw(a::Number,b::Number,c::Number)
	v(z) = SVector(1,real(z),imag(z))
	det([v(a) v(b) v(c)]) > 0
end

function fdtangent(z,t::Real) 
	ϵ = eps(typeof(float(t)))
	ϵ3 = ϵ^(1/3)
	# Use a finite difference approximation: 1st order at edges, 2nd otherwise
	if t < ϵ3
		t0, t1 = t, t+sqrt(ϵ)
	elseif t > 1-ϵ3
		t0, t1 = t-sqrt(ϵ), t
	else
		t0, t1 = t-ϵ3,t+ϵ3 
	end
	return (z(t1)-z(t0))/(t1-t0)
end

function adaptpoints(point,utangent,a,b;depth=6,curvemax=0.05)
	function refine(tl,tr,zl,zr,τl,τr,maxdz,d=depth)
		# approximately the stepsize over radius of curvature
		dzkap = dist(τr,τl)

		tm = (tl+tr)/2
		zm = point(tm)
		τm = utangent(tm) 

		if d > 0 && (dzkap > curvemax || dist(zr,zl) > maxdz )
			zl = refine(tl,tm,zl,zm,τl,τm,maxdz,d-1)
			zr = refine(tm,tr,zm,zr,τm,τr,maxdz,d-1)
			return [zl;zm;zr]
		else
			return zm
		end
	end


	d = (b-a)/4
	tt = d*[0,0.196,0.41,0.592,0.806] 
	t = [a .+ tt; a + d .+ tt; a + 2d .+ tt; a + 3d .+ tt; b]
	z = point.(t)
	τ = utangent.(t)

	if z[1] isa Spherical
		dist = (u,v) -> norm(S2coord(u)-S2coord(v))
	else
		dist = (u,v) -> abs(u-v)
	end
	m = length(t) 
	scale = maximum(dist(z[i],z[j]) for i=2:m-1, j=2:m-1 if j > i)
	zfinal = z[[1]]
	for j = 1:length(t)-1
		znew = refine(t[j],t[j+1],z[j],z[j+1],τ[j],τ[j+1],scale/25)
		append!(zfinal,znew)
		push!(zfinal,z[j+1])
	end
	return zfinal
end