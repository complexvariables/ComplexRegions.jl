# Type and constructors
"""
	(type) Circle{T<:AnyComplex} in the complex plane

Each `Circle` type is parameterized according to the common type of its complex input arguments.
"""
struct Circle{T} <: AbstractClosedCurve{T}
    center::AnyComplex{T}
    radius::T
    ccw::Bool
	function Circle{T}(z, r, ccw=true) where {T}
		Tz = complex(convert_real_type(T, z))
		Tr = convert(T, r)
		return new(Tz, Tr, ccw)
	end
end
"""
	Circle(zc, r, ccw=true)

Construct the circle with given center `zc`, radius `r`, and orientation (defaults to counterclockwise).

	Circle(a, b, c)

Construct the circle passing through the given three numbers. Orientation is determined so that the values are visited in the given order. If the three points are collinear (including when one of the given values is infinite), a `Line` is returned instead.
"""
Circle(z::AnyComplex, r::Real, ccw::Bool=true) = Circle{real_type(z)}(z, r, ccw)
Circle(z::Number, r::Real, ccw::Bool=true) = Circle(complex(float(z)), r, ccw)

# Construction by three points
function Circle(a::Number, b::Number, c::Number)
    a, b, c = promote(complex.(float.((a, b, c)))...)
	T = real_type(a)
    isinf(a) && return Line(b, c)
    isinf(b) && return Line(c, a)
    isinf(c) && return Line(a, b)
    # Use intersection of chord bisectors to find the center of the circle.
    w = (a - c) / 2
    d1, d2 = a - b, c - b
    M = SMatrix{2,2}(real(d1), imag(d1), real(d2), imag(d2))
    if abs(det(M)) < tolerance(T)
        # Essentially collinear points
        return Line(a, b)
    else
        p = M \ SVector(imag(w), -real(w))
        cen = (a + b) / 2 - 1im * p[1] * d1
        ccw = isccw(a - cen, b - cen, c - cen)
        return Circle{T}(cen, abs(a - cen), ccw)
    end
end

function Base.convert(::Type{Circle{T}}, C::Circle{S}) where {T,S}
	return Circle{T}(convert_real_type(T, C.center), convert(T, C.radius), C.ccw)
end
convert_real_type(::Type{T}, C::Circle{S}) where {T<:Real,S} = convert(Circle{T}, C)
Base.promote_rule(::Type{<:Circle{T}}, ::Type{<:Circle{S}}) where {T,S} = Circle{promote_type(T,S)}

# Required methods
function point(C::Circle, t::Real)
    z = C.ccw ? cispi(2t) : cispi(-2t)
    C.center + typeof(C.center)(z * C.radius)
end

ispositive(C::Circle) = C.ccw
arclength(C::Circle{T}) where T = 2T(π) * C.radius

"""
	arg(C::Circle,z)

Find the parameter argument `t` such that `C(t)==z` is true.

This gives undefined results if `z` is not actually on the circle.
"""
function arg(C::Circle{T}, z::Number) where T
    α = angle(z - C.center) / (2T(π))
    C.ccw ? mod(α, 1) : mod(-α, 1)
end

function unittangent(C::Circle, t::Real)
    C.ccw ? 1im * cispi(2t) : -1im * cispi(-2t)
end
function tangent(C::Circle, t::Real)
    2*(π * C.radius * unittangent(C, t))
end

Base.:+(C::Circle, z::Number) = Circle(C.center + z, C.radius, C.ccw)
Base.:-(C::Circle) = Circle(-C.center, C.radius, C.ccw)
Base.:*(C::Circle, z::Number) = Circle(C.center * z, C.radius * abs(z), C.ccw)

"""
	inv(C)
Invert the circle `C` through the origin. In general the inverse is a `Circle`, though the result is a `Line` if `C` passes through the origin.
"""
function inv(C::Circle{T}) where T
    w = 1 ./ point(C, SVector(T(0), T(1)/4, T(1)/2))
    Circle(w...)
end

# Other methods
isfinite(::Circle) = true
conj(C::Circle) = Circle(conj(C.center), C.radius, !C.ccw)
reverse(C::Circle) = Circle(C.center, C.radius, !C.ccw)

"""
	isapprox(C1::Circle,C2::Circle; tol=<default>)
	C1 ≈ C2
Determine if `C1` and `C2` represent the same circle, irrespective of the type or values of its parameters. Identity is determined by agreement within `tol`, which is interpreted as the weaker of absolute and relative differences.
"""
function isapprox(C1::Circle{S}, C2::Circle{T}; tol=tolerance(S, T)) where {S,T}
    return isapprox(C1.center, C2.center, rtol=tol, atol=tol) &&
           isapprox(C1.radius, C2.radius, rtol=tol, atol=tol)
end

function winding(C::Circle, z::Number)
    w = abs(z - C.center) < C.radius ? 1 : 0
    C.ccw ? w : -w
end

"""
	dist(z,C::Circle)
Compute the distance from number `z` to the circle `C`.
"""
dist(z::Number, C::Circle) = abs(abs(z - C.center) - C.radius)

"""
	closest(z,C::Circle)
Find the point on circle `C` that lies closest to `z`.
"""
closest(z::Number, C::Circle) = C.center + C.radius * sign(z - C.center)

"""
	reflect(z,C::Circle)
Reflect the value `z` across the circle `C`. (For reflection of a circle through a point, use translation and negation.)
"""
function reflect(z::Number, C::Circle)
    ζ = z - C.center
    ζ == 0 ? convert(typeof(float(z)), Inf) : C.center + C.radius^2 / conj(ζ)
end

unitcircle = Circle(0, 1)

# COV_EXCL_START
function show(io::IO, C::Circle)
    orient = C.ccw ? "ccw" : "cw"
    print(IOContext(io, :compact => true), "Circle(", C.center, ",", C.radius, ",", orient, ")")
end

function show(io::IO, ::MIME"text/plain", C::Circle{T}) where {T}
    orient = C.ccw ? "positively" : "negatively"
    print(io, "Circle{$T} in the complex plane:\n   centered at (", C.center, ") with radius $(C.radius), $(orient) oriented")
end
# COV_EXCL_END
