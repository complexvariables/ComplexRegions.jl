const LineOrCircle{T} = Union{Line{T},Circle{T}}

"""
	(type) Möbius
Representation of a Möbius or bilinear transformation.
"""
struct Möbius <: Function
    # interpreted as [a,b,c,d], where f(z)=(az+b)/(cz+d)
    coeff::SVector{4}
end

#
# Construction
#

const Mobius = Möbius   # for us lazy Americans
"""
	Möbius(a, b, c, d)
Construct the `Möbius` map ``z ↦ (az+b) / (cz+d)`` from its coefficients.
"""
Möbius(a::Number, b::Number, c::Number, d::Number) = Möbius(SVector(a, b, c, d))

"""
	Möbius(A::AbstractMatrix)
Construct the `Möbius` map ``z ↦ (az+b) / (cz+d)`` from the matrix `[a b;c d]`.
"""
Möbius(A::AbstractMatrix) = Möbius(SVector(A[1, 1], A[1, 2], A[2, 1], A[2, 2]))

"""
	Möbius(z::AbstractVector, w::AbstractVector)
Construct the `Möbius` map that transforms the points `z[k]` to `w[k]` for k=1,2,3.
Values of `Inf` are permitted in both vectors.
"""
function Möbius(source::AbstractVector, image::AbstractVector)
    # finds the coeffs of map from (0, 1, Inf) to given points
    function standard(x, y, z)
        if isinf(x)
            return z, y - z, 1, 0
        elseif isinf(y)
            return -z, x, -1, 1
        elseif isinf(z)
            return y - x, x, 0, 1
        else
            xy, yz = y - x, z - y
            return z * xy, x * yz, xy, yz
        end
    end

    za, zb, zc, zd = standard(source...)
    wa, wb, wc, wd = standard(image...)
    A = SMatrix{2,2}(wa, wc, wb, wd) * SMatrix{2,2}(zd, -zc, -zb, za)
    return Möbius(A[1, 1], A[1, 2], A[2, 1], A[2, 2])
end

"""
	Möbius(C1, C2)
Construct a `Möbius` map that transforms the curve `C1` to `C2`. Both curves must be either a `Line` or `Circle`. (These maps are not uniquely determined.)
"""
function Möbius(c1::LineOrCircle{T}, c2::LineOrCircle{T}) where T
	t = SVector(T(1), T(1)/4, T(1)/2)
	return Möbius(points(c1, t), points(c2, t))
end

#
# Evaluation
#

"""
	f(z::Number)
Evaluate the `Möbius` map `f` at a real or complex value `z`.
"""
function (f::Möbius)(z::Number)
    a, b, c, d = f.coeff
    if isinf(z)
        num, den = a, c
    else
        num, den = a * z + b, c * z + d
    end
    # note: 1 over complex zero is NaN, not Inf
    iszero(den) ? complex(abs(num) / abs(den)) : num / den
end

"""
	f(C::Union{Circle,Line})
Find the image of the circle or line `C` under the `Möbius` map `f`. The result is also either a `Circle` or a `Line`.
"""
(f::Möbius)(C::LineOrCircle{T}) where T = Circle(f.(points(C, [T(0), T(1)/4, T(1)/2]))...)

"""
	f(C::Union{Arc,Segment})
Find the image of the arc or segment `C` under the `Möbius` map `f`. The result is also either an `Arc` or a `Segment`.
"""
(f::Möbius)(C::Union{Arc{T}, Segment{T}}) where T = Arc(f.(points(C, [T(0), T(1)/2, T(1)]))...)

"""
	f(R::Union{Disk,Halfplane})
If `R` is an `Disk` or an `Halfplane`, find its image under the `Möbius` map `f`. The result is also either an `Disk` or an `Halfplane`.
# Examples
```julia-repl
julia> f = Möbius(Line(-1,1),Circle(0,1))
Möbius transformation:

   (1.0 + 0.9999999999999999im) z + (3.666666666666666 - 1.666666666666667im)
   ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
   (1.0 + 0.9999999999999999im) z + (-1.666666666666666 + 3.6666666666666665im)

julia> f(upperhalfplane)
Disk interior to:
   Circle(-5.55112e-17+2.22045e-16im,1.0)

julia> isapprox(ans,unitdisk)
true
```
"""
function (f::Möbius)(R::Union{
							Disk,
							ExteriorSimplyConnectedRegion,
							Halfplane}
					)
	return interior(f(R.boundary))
end

"""
	inv(f::Möbius)
Find the inverse of a `Möbius` transformation. This is the functional inverse, not 1/f(z).
"""
inv(f::Möbius) = Möbius(f.coeff[4], -f.coeff[2], -f.coeff[3], f.coeff[1])

"""
	∘(f::Möbius,g::Möbius)
Compose two `Möbius` transformations.
"""
function ∘(f::Möbius, g::Möbius)
    A = SMatrix{2,2}(f.coeff[1], f.coeff[3], f.coeff[2], f.coeff[4])
    B = SMatrix{2,2}(g.coeff[1], g.coeff[3], g.coeff[2], g.coeff[4])
    C = A * B
    Möbius(C[1, 1], C[1, 2], C[2, 1], C[2, 2])
end

#
# Display
#

# COV_EXCL_START
function Base.show(io::IO, f::Möbius)
    a, b, c, d = f.coeff
    print(IOContext(io, :compact => true), "Möbius map z --> (($a)*z + $b) / (($c)*z + $d)")
end

function Base.show(io::IO, ::MIME"text/plain", f::Möbius)
    a, b, c, d = [repr("text/plain", x) for x in f.coeff]
    numer = "(" * a * ") z + (" * b * ")"
    denom = "(" * c * ") z + (" * d * ")"
    hline = reduce(*, fill("–", max(length(numer), length(denom))))
    print(io, "Möbius transformation:\n\n   ", numer, "\n   ", hline, "\n   ", denom)
end
# COV_EXCL_STOP
