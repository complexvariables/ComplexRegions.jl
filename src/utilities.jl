# Scale from/to [0,1]. 
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

# Are three given points arranged in counterclockwise order? 
function isccw(a::Number,b::Number,c::Number)
	v(z) = SVector(1,real(z),imag(z))
	det([v(a) v(b) v(c)]) > 0
end

# Use 2nd order finite differences to approximate a tangent.
function fdtangent(z,t::Real) 
	ϵ = eps(typeof(float(t)))
	ϵ3 = 0.5*ϵ^(1/3)
	if t < ϵ3
		τ = (-1.5*z(t) + 2*z(t+ϵ3) - 0.5*z(t+2ϵ3)) / ϵ3
	elseif t > 1-ϵ3
		τ = (1.5*z(t) - 2*z(t-ϵ3) + 0.5*z(t-2ϵ3)) / ϵ3
	else
		τ = (z(t+ϵ3) - z(t-ϵ3)) / (2ϵ3)
	end
	return τ
end

# Select points adaptively to make a smooth-appearing curve. 
function adaptpoints(point,utangent,a,b;depth=6,curvemax=0.05)
	function refine(tl,tr,zl,zr,τl,τr,maxdz,d=depth)		
		dzkap = dist(τr,τl)  # approximately, the stepsize over radius of curvature
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
	tt = d*[0,0.196,0.41,0.592,0.806]   # avoid common symmetry points
	t = [a .+ tt; a + d .+ tt; a + 2d .+ tt; a + 3d .+ tt; b]
	z = point.(t)
	τ = utangent.(t)

	# on the Riemann sphere, use distance in R^3
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

# Do adaptive integration to estimate the integral of `f` over [`a`,`b`] to desired
# error tolerance `tol`.
function intadapt(f,a,b,tol)
    # Use error estimation and recursive bisection.
    function do_integral(a,fa,b,fb,m,fm,tol,depth)
        # These are the two new nodes and their f-values.
        xl = (a+m)/2;  fl = f(xl);
        xr = (m+b)/2;  fr = f(xr);
        t = [a,xl,m,xr,b]              # all 5 nodes at this level

        # Compute the trapezoid values iteratively.
        h = (b-a)
        T = [0.,0.,0.]
        T[1] = h*(fa+fb)/2
        T[2] = T[1]/2 + (h/2)*fm
        T[3] = T[2]/2 + (h/4)*(fl+fr)

        S = (4*T[2:3]-T[1:2]) / 3      # Simpson values
        E = (S[2]-S[1]) / 15           # error estimate

        if abs(E) < tol*(1+abs(S[2]))  # acceptable error?
            Q = S[2]                   # yes--done
		else
			if depth==0
				@warn "Too many recursions to determine integral"
				Q = S[2] 
			else
    	        # Error is too large--bisect and recurse.
        	    QL = do_integral(a,fa,m,fm,xl,fl,tol,depth-1)
            	QR = do_integral(m,fm,b,fb,xr,fr,tol,depth-1)
				Q = QL + QR
			end
        end
        return Q
    end

    m = (b+a)/2
    Q = do_integral(a,f(a),b,f(b),m,f(m),tol,50)
    return Q
end

function enclosing_circle(z::AbstractVector,expansion=2)
	xa,xb = extrema(real(z))
	ya,yb = extrema(imag(z))
	zc = complex((xa+xb)/2,(ya+yb)/2)
	R = length(z) > 1 ? maximum(@. abs(z - zc)) : max(1,abs(zc))
	return zc,expansion*R
end

function enclosing_box(z::AbstractVector,expansion=2)
	zc = sum(z)/length(z)
	dz = z .- zc
	rx = length(z) > 1 ? maximum(@. abs(real(dz))) : max(1,abs(real(zc)))
	ry = length(z) > 1 ? maximum(@. abs(imag(dz))) : max(1,abs(imag(zc)))
	return real(zc).+expansion*[-rx,rx],imag(zc).+expansion*[-ry,ry] 
end

# indices of the closest pair of points from two lists
function argclosest(z1,z2)
	i1 = [argmin(abs.(z1.-z)) for z in z2]
	i2 = argmin(abs.(z2.-z1[i1]))
	return i1[i2],i2
end