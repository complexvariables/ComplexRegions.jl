abstract type AbstractMap end

struct Möbius <: AbstractMap 
	coeff::SVector{4} 
end 

const Mobius = Möbius 
Möbius(a::Number,b::Number,c::Number,d::Number) = Möbius(SVector(a,b,c,d))
Line_Circle = Union{Line,Circle}

function Möbius(source::AbstractVector,image::AbstractVector) 
	# finds the coeffs of map from (0,1,Inf) to given points 
	function standard(x,y,z) 
		if isinf(x)
			return z,y-z,1,0
		elseif isinf(y)
			return -z,x,-1,1
		elseif isinf(z) 
			return y-x,x,0,1
		else
			xy,yz = y-x,z-y
			return z*xy,x*yz,xy,yz
		end
	end

	za,zb,zc,zd = standard(source...)
	wa,wb,wc,wd = standard(image...)
	A = SMatrix{2,2}(wa,wc,wb,wd)*SMatrix{2,2}(zd,-zc,-zb,za)
	return Möbius(A[1,1],A[1,2],A[2,1],A[2,2])
end

Möbius(c1::Line_Circle,c2::Line_Circle) = Möbius(point(c1,[0,0.25,0.5]),point(c2,[0,0.25,0.5]))

#
# Evaluation
#
# at point
function (f::Möbius)(z::Number)
	a,b,c,d = f.coeff
	if isinf(z) 
		num,den = a,c
	else 
		num,den = a*z+b,c*z+d
	end
	# note: 1 over complex zero is NaN, not Inf
	iszero(den) ? complex(abs(num)/abs(den)) : num/den
end

# for Circle or Line 
(f::Möbius)(C::Line_Circle) = Circle( f.(point(C,[0,0.25,0.5]))... )
(f::Möbius)(C::Union{Arc,Segment}) = Arc( f.(point(C,[0,0.5,1]))... )

# for Disk or Halfplane
(f::Möbius)(R::Union{AbstractDisk,AbstractHalfplane}) = interior(f(R.boundary)) 

inv(f::Möbius) = Möbius(f.coeff[4],-f.coeff[2],-f.coeff[3],f.coeff[1])

