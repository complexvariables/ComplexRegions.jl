# Type  
"""
	(type) Ray{T<:AnyComplex} in the complex plane 
Each `Ray` type is parameterized according to the common type of its input arguments. 
"""
struct Ray{T<:AnyComplex} <: AbstractCurve 
	base::T 
	angle::AbstractFloat  
	reverse::Bool
	function Ray{T}(a,d,rev=false) where T<:AnyComplex
		new(a,mod2pi(d),rev)
	end
end

# Untyped constructor
"""
	Ray(a,z,reverse=false)
Construct the ray starting at `a` and extending to infinity in the direction parallel to `z`. If `reverse` is true, the ray is considered to extend from infinity to `a` along `-z`.
"""
function Ray(a::Number,d::Number,rev=false) 
	a = complex(float(a))
	Ray{typeof(a)}(a,float(d),rev)
end

# Complex type converters
for ctype in [:Spherical,:Polar,:Complex]
	docstr = """
	$ctype(::Ray)
Convert to `Ray{$ctype}`. This is useful for plotting curves in a desired way.
"""
	@eval begin 
		@doc $docstr ->
		function $ctype(R::Ray{T}) where T<:AnyComplex 
			Ray($ctype(R.base),R.angle,R.reverse)
		end
	end
end	

# Required methods
arclength(R::Ray) = Inf
function point(R::Ray{T},t::Real) where T
	if R.reverse 
		t = 1-t 
	end
	# avoid NaNs 
	if t==0 
		R.base 
	elseif t==1 
		T(Inf)
	else
		R.base + t/(1-t)*exp(complex(0,R.angle))
	end	
end
(C::Ray)(t::Real) = point(C,t)
""" 
	arg(R::Ray,z) 
Find the parameter argument `t` such that `R(t)==z` is true. 

This gives undefined results if `z` is not actually on the ray. 
"""
function arg(R::Ray,z::Number)
	if isinf(z)
		t = 1
	else
		δ = abs(z - R.base) 
		t = δ / (1+δ)
	end
	return R.reverse ? 1-t : t 
end
tangent(R::Ray,t::Real) = tangent(R::Ray)
tangent(R::Ray) = R.reverse ? -exp(1im*R.angle) : exp(1im*R.angle)

# Other methods
isfinite(::Ray) = false
conj(R::Ray) = Ray(conj(R.base),-R.angle,R.reverse)
reverse(R::Ray) = Ray(R.base,R.angle,!R.reverse)
"""
	R + z
	z + R 
Translate the ray `R` by a number `z`. 
"""
+(R::Ray,z::Number) = Ray(R.base+z,R.angle,R.reverse)
+(z::Number,R::Ray) = Ray(R.base+z,R.angle,R.reverse)
"""
	R - z
Translate the ray `R` by a number `-z`.

	-R 
	z - R 
Negate a ray `R` (reflect through the origin), and optionally translate by a number `z`.
"""
-(R::Ray,z::Number) = Ray(R.base-z,R.angle,R.reverse)
-(z::Number,R::Ray) = Ray(z-R.base,R.angle,R.reverse)
-(R::Ray) = Ray(-R.base,mod2pi(R.angle+pi),R.reverse)
# these need to recompute the final parameter values
"""
	z*R 
	R*z 
Multiply the ray `R` by real or complex number `z`; i.e., scale and rotate it about the origin.
"""
*(R::Ray,z::Number) = Ray(z*R.base,mod2pi(R.angle+sign(z)),R.reverse)
*(z::Number,R::Ray) = *(R,z)
"""
	R/z 
Multiply the ray `R` by the number `1/z`; i.e., scale and rotate it about the origin.

	z/R 
	inv(R) 
Invert the ray `R` through the origin (and optionally multiply by the number `1/z`). In general the inverse is an `Arc`.
"""
/(R::Ray,z::Number) = *(R,1/z)
function /(z::Number,R::Ray) 
	w = z./point(R,[0,0.5,1])
	Arc(w...)
end
inv(R::Ray) = 1/R

"""
	isapprox(R1::Ray,R2::Ray; tol=<default>) 
	R1 ≈ R2 
Determine if `R1` and `R2` represent the same segment, irrespective of the type or values of its parameters. Identity is determined by agreement within `tol`, which is interpreted as the weaker of absolute and relative differences.
"""
function isapprox(R1::Ray,R2::Ray;tol=DEFAULT[:tol])
	return isapprox(R1.base,R2.base,rtol=tol,atol=tol) && (abs(mod2pi(R1.angle-R2.angle)) < tol)
end

""" 
	isleft(z,R::Ray) 
Determine whether the number `z` lies "to the left" of ray `R`. This means that the angle it makes with `tangent(R)` is in the interval (0,π).

Note that `isleft` and `isright` are *not* logical opposites; a point on the (extended) ray should give `false` in both cases.
"""
function isleft(z::Number,R::Ray) 
	a,b = point(R,[0.2,0.8])  # accounts for reversal
	(real(b)-real(a)) * (imag(z)-imag(a)) > (real(z)-real(a)) * (imag(b)-imag(a))
end
""" 
	isright(z,R::Ray) 
Determine whether the number `z` lies "to the right" of ray `R`. This means that the angle it makes with `tangent(R)` is in the interval (-π,0).

Note that `isleft` and `isright` are *not* logical opposites; a point on the (extended) ray should give `false` in both cases.
"""
function isright(z::Number,R::Ray) 
	a,b = point(R,[0.2,0.8])  # accounts for reversal
	(real(b)-real(a)) * (imag(z)-imag(a)) < (real(z)-real(a)) * (imag(b)-imag(a))
end

""" 
	dist(z,R::Ray) 
Compute the distance from number `z` to the ray `R`. 
"""
dist(z::Number,R::Ray) = abs(z - closest(z,R))
""" 
	closest(z,R::Ray) 

Find the point on ray `R` that lies closest to `z`.
"""
function closest(z::Number,R::Ray) 
	# translate and rotate to positive Re axis
	s = exp(complex(0,R.angle))
	ζ = (z-R.base)/s
	R.base + max(real(ζ),0)*s
end

sign(R::Ray) = tangent(R)

# Display methods
function show(io::IO,R::Ray{T}) where {T}
	print(IOContext(io,:compact=>true),"Ray(",R.base,",",R.angle,",",R.reverse,")")
end

function show(io::IO,::MIME"text/plain",R::Ray{T}) where {T}
	if R.reverse 
		print(io,"Ray to ",R.base," at angle ",R.angle)
	else
		print(io,"Ray from ",R.base," at angle ",R.angle)
	end
end

# in the plane, use two (finite) points for plotting
plotdata(R::Ray{T}) where T<:Union{Complex,Polar} = R.reverse ? [R(0.3),R.base] : [R.base,R(0.7)]
