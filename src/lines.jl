# Type and constructors
# Default constructor is for the line through two points
"""
	(type) Line{T<:AnyComplex} in the complex plane

Each `Line` type is parameterized according to the common type of its complex input arguments.
"""
struct Line{T} <: AbstractClosedCurve{T}
    base::Union{T,AnyComplex{T}}
    direction::Union{T,AnyComplex{T}}
    Line{T}(a, b) where {T} = new(convert_real_type(T, a), sign(convert_real_type(T, b - a)))
end
"""
	Line(a, b)

Construct the line passing through the two given points.

	Line(a, direction=z)

Construct the line through the given point and parallel to the complex sign of the given `direction` value. """
function Line(a::Number, b::Number)
    a, b = promote(complex(float(a)), complex(float(b)))
    Line{real_type(a)}(a, b)
end

# Construct by giving a keyword
Line(z::Number; direction::Number) = Line(z, z + direction)

# # Complex type converters
# for ctype in [:Spherical, :Polar, :Complex]
#     @eval begin
#         function $ctype(L::Line{T}) where {T<:AnyComplex}
#             Line($ctype(L.base), $ctype(L.base + L.direction))
#         end
#     end
# end

function Base.convert(::Type{Line{T}}, L::Line{S}) where {T,S}
	return Line{T}(convert_real_type(T, L.base), convert_real_type(T, L.direction))
end
convert_real_type(::Type{T}, L::Line{S}) where {T<:Real,S} = convert(Line{T}, L)
Base.promote_rule(::Type{<:Line{T}}, ::Type{<:Line{S}}) where {T,S} = Line{promote_type(T,S)}

# Required methods
arclength(::Line) = Inf
point(L::Line, t::Real) = L.base + (2t - 1) / (t - t^2) * L.direction

ispositive(C::Line) = true  # this seems arbitrary...

"""
	arg(L::Line, z)

Find a parameter argument `t` such that `L(t)==z` is true. For an infinite `z`, return zero (but note that `L(1)` is also infinity).

This gives undefined results if `z` is not actually on the line.
"""
function arg(L::Line, z::Number)
    isinf(z) && return zero(z)
    del = z - L.base
    if abs(real(del)) > abs(imag(del))
        del = real(del)
        α = real(L.direction)
    else
        del = imag(del)
        α = imag(L.direction)
    end
    for t in realroots(del, 2α - del, -α)
        0 - eps(del) ≤ t ≤ 1 + eps(del) && return t
    end
    return nothing
end

unittangent(L::Line, t::Real=0) = L.direction
tangent(L::Line, t::Real) = (1 / t^2 + 1 / (1 - t)^2) * L.direction

Base.:+(L::Line, z::Number) = Line(L.base + z, direction=L.direction)
Base.:-(L::Line) = Line(-L.base, direction=-L.direction)
Base.:*(L::Line, z::Number) = Line(L.base * z, direction=L.direction * sign(z))

"""
	inv(L)
Invert a line `L` through the origin. In general the inverse is a `Circle` through the inverse of any three points on the line.
"""
function inv(L::Line{T}) where T
    w = 1 ./ point(L, SVector(T(1) / 4, T(1) / 2, T(3) / 4))
    Circle(w...)
end

# Other methods
isfinite(::Line) = false
slope(L::Line) = imag(L.direction) / real(L.direction)
angle(L::Line) = angle(L.direction)
conj(L::Line) = Line(conj(L.base), direction=conj(L.direction))
reverse(L::Line) = Line(L.base, direction=-L.direction)

"""
	isapprox(L1::Line, L2::Line; tol=<default>)
	L1 ≈ L2

Determine if `L1` and `L2` represent the same line, irrespective of the type or values of its parameters. Identity is determined by agreement within `tol`, which is interpreted as the weaker of absolute and relative differences.
"""
function isapprox(L1::Line{T}, L2::Line{S}; tol=tolerance(S, T)) where {T,S}
    dz = L1.base - L2.base
    w1, w2 = L1.direction, L2.direction
    return isapprox(real(w1) * imag(w2), imag(w1) * real(w2), rtol=tol, atol=tol) &&
           isapprox(real(w1) * imag(dz), imag(w1) * real(dz), rtol=tol, atol=tol)
end

"""
	isleft(z, L::Line)

Determine whether the number `z` lies "to the left" of line `L`. This means that the angle it makes with `tangent(L)` is in the interval (0,π).

Note that `isleft` and `isright` are *not* logical opposites; a point on the curve should give `false` in both cases.
"""
isleft(z::Number, L::Line) = π > angle((z - L.base) / L.direction) > 0

"""
	isright(z, L::Line)

Determine whether the number `z` lies "to the right" of line `L`. This means that the angle it makes with `tangent(L)` is in the interval (-π,0).

Note that `isleft` and `isright` are *not* logical opposites; a point on the curve should give `false` in both cases.
"""
isright(z::Number, L::Line) = -π < angle((z - L.base) / L.direction) < 0

winding(L::Line, z::Number) = isleft(z, L) ? 1 : 0

"""
	dist(z, L::Line)

Compute the distance from number `z` to the line `L`.
"""
dist(z::Number, L::Line) = imag((z - L.base) / L.direction)

"""
	closest(z, L::Line)

Find the point on line `L` that lies closest to `z`.
"""
function closest(z::Number, L::Line)
    s = L.direction
    L.base + real((z - L.base) / s) * s
end

"""
	reflect(z, L::Line)

Reflect the value `z` across the line `L`. (For reflection of a line through a point, use translation and negation.)
"""
function reflect(z::Number, L::Line)
    ζ = z - L.base
    L.base + L.direction * conj(ζ / L.direction)
end

function Base.show(io::IO, L::Line)
    print(IOContext(io, :compact => true), "Line(...", L(0.5), "...", L((sqrt(5) - 1) / 2), "...)")
end

function Base.show(io::IO, ::MIME"text/plain", L::Line{T}) where {T}
    print(io, "Line{$T} in the complex plane:\n   through (", L.base, ") parallel to (", L.direction, ")")
end

# Two points are enough to draw a line (though not on the sphere), and we want to avoid infinity.
plotdata(L::Line) = point(L, [0.1, 0.9])
