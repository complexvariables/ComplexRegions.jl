# Type  
struct Ray{T<:AnyComplex} <: AbstractCurve 
	base::T 
	angle::AbstractFloat  
	reverse::Bool
	function Ray{T}(a,d,rev=false) where T<:AnyComplex
		new(a,mod2pi(d),rev)
	end
end

# Untyped constructor
function Ray(a::Number,d::Number,rev=false) 
	a = complex(float(a))
	Ray{typeof(a)}(a,float(d),rev)
end

# Complex type converters
for ctype in [:Spherical,:Polar,:Complex]
	@eval begin 
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
isbounded(::Ray) = false
conj(R::Ray) = Ray(conj(R.base),-R.angle,R.reverse)
reverse(R::Ray) = Ray(R.base,R.angle,!R.reverse)
+(R::Ray,z::Number) = Ray(R.base+z,R.angle,R.reverse)
+(z::Number,R::Ray) = Ray(R.base+z,R.angle,R.reverse)
-(R::Ray,z::Number) = Ray(R.base-z,R.angle,R.reverse)
-(z::Number,R::Ray) = Ray(z-R.base,R.angle,R.reverse)
-(R::Ray) = Ray(-R.base,mod2pi(R.angle+pi),R.reverse)
# these need to recompute the final parameter values
*(R::Ray,z::Number) = Ray(z*R.base,mod2pi(R.angle+sign(z)),R.reverse)
*(z::Number,R::Ray) = *(R,z)
/(R::Ray,z::Number) = *(R,1/z)
function /(z::Number,R::Ray) 
	w = z./point(R,[0,0.5,1])
	Arc(w...)
end
inv(R::Ray) = 1/R

function isapprox(R1::Ray,R2::Ray;tol=1e-12)
	return isapprox(R1.base,R2.base,rtol=tol,atol=tol) && (abs(mod2pi(R1.angle-R2.angle)) < tol)
end

dist(z::Number,R::Ray) = abs(z - closest(z,R))
function closest(z::Number,R::Ray) 
	# translate and rotate to positive Re axis
	s = exp(complex(0,R.angle))
	ζ = (z-R.base)/s
	R.base + max(real(ζ),0)*s
end

sign(R::Ray) = tangent(R)

function isleft(z::Number,R::Ray) 
	a,b = point(R,[0.2,0.8])  # accounts for reversal
	(real(b)-real(a)) * (imag(z)-imag(a)) > (real(z)-real(a)) * (imag(b)-imag(a))
end

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

plotdata(R::Ray{T}) where T<:Union{Complex,Polar} = R.reverse ? [R(0.3),R.base] : [R.base,R(0.7)]
