# Type
"""
	(type) Segment{T<:AbstractFloat} in the complex plane

Each `Segment` type is parameterized according to the common type of its complex input arguments.
"""
struct Segment{T} <: AbstractCurve{T}
    za::Union{T,AnyComplex{T}}
    zb::Union{T,AnyComplex{T}}
    function Segment{T}(a, b) where {T}
        @assert isfinite(a) && isfinite(b)
        new(convert_real_type(T, a), convert_real_type(T, b))
    end
end

# Untyped constructor
"""
	Segment(a, b)

Consruct a segment that starts at value `a` and ends at `b`.
"""
function Segment(a::Number, b::Number)
	T = promote_type(real_type(float(a)), real_type(float(b)))
	a, b = convert_real_type(T, a), convert_real_type(T, b)
	Segment{T}(a, b)
end

# We want a real segment to have real values
function Segment(a::T, b::S) where {T<:Real,S<:Real}
	a, b = promote(float(a), b)
	Segment{typeof(a)}(a, b)
end

function Base.convert(::Type{Segment{T}}, L::Segment{S}) where {T,S}
	return Segment{T}(convert_real_type(T, L.za), convert_real_type(T, L.zb))
end
convert_real_type(::Type{T}, L::Segment{S}) where {T<:Real,S} = convert(Segment{T}, L)
Base.promote_rule(::Type{<:Segment{T}}, ::Type{<:Segment{S}}) where {T,S} = Segment{promote_type(T,S)}

# Required methods
arclength(S::Segment) = abs(S.zb - S.za)
point(S::Segment, t::Real) = (1 - t) * S.za + t * S.zb

"""
	arg(S::Segment,z)

Find the parameter argument `t` such that `S(t)==z` is true.

This gives undefined results if `z` is not actually on the segment.
"""
arg(S::Segment, z::Number) = real((z - S.za) / (S.zb - S.za))
tangent(S::Segment) = S.zb - S.za
tangent(S::Segment, t::Real) = tangent(S)
unittangent(S::Segment) = sign(S.zb - S.za)
unittangent(S::Segment, t::Real) = unittangent(S)
normal(S::Segment) = normal(S, 0)

Base.:+(S::Segment, z::Number) = Segment(S.za + z, S.zb + z)
Base.:-(S::Segment) = Segment(-S.za, -S.zb)
Base.:*(S::Segment, z::Number) = Segment(S.za * z, S.zb * z)

"""
	inv(S)
Invert the segment `S` through the origin. In general the inverse is an `Arc`, though the result is a `Segment` if `S` would pass through the origin when extended.
"""
function inv(S::Segment)
    w = 1 ./ point(S, SVector(0, 1//2, 1))
    return Arc(w...)
end

# Other methods
isfinite(::Segment) = true
conj(S::Segment) = Segment(conj(S.za), conj(S.zb))
reverse(S::Segment) = Segment(S.zb, S.za)
Base.isreal(S::Segment) = isreal(S.za) && isreal(S.za)
sign(S::Segment) = unittangent(S)

"""
	isapprox(S1::Segment,S2::Segment; tol=<default>)
	S1 ≈ S2

Determine if `S1` and `S2` represent the same segment, irrespective of the type or values of its parameters. Identity is determined by agreement within `tol`, which is interpreted as the weaker of absolute and relative differences.
"""
function isapprox(S1::Segment{S}, S2::Segment{T}; tol=tolerance(S, T)) where {S,T}
    return isapprox(S1.za, S2.za, rtol=tol, atol=tol) &&
           isapprox(S1.zb, S2.zb, rtol=tol, atol=tol)
end

"""
	isleft(z,S::Segment)

Determine whether the number `z` lies "to the left" of segment `S`. This means that the angle it makes with `tangent(S)` is in the interval (0,π).

Note that `isleft` and `isright` are *not* logical opposites; a point on the (extended) segment should give `false` in both cases.
"""
function isleft(z::Number, S::Segment)
    a, b = S.za, S.zb
    return (real(b) - real(a)) * (imag(z) - imag(a)) > (real(z) - real(a)) * (imag(b) - imag(a))
end
"""
	isright(z,S::Segment)

Determine whether the number `z` lies "to the right" of segment `S`. This means that the angle it makes with `tangent(S)` is in the interval (-π,0).

Note that `isleft` and `isright` are *not* logical opposites; a point on the (extended) segment should give `false` in both cases.
"""
function isright(z::Number, S::Segment)
    a, b = S.za, S.zb
    return (real(b) - real(a)) * (imag(z) - imag(a)) < (real(z) - real(a)) * (imag(b) - imag(a))
end

"""
	dist(z,S::Segment)

Compute the distance from number `z` to the segment `S`.
"""
dist(z::Number, S::Segment) = abs(z - closest(z, S))

"""
	closest(z,S::Segment)

Find the point on segment `S` that lies closest to `z`.
"""
function closest(z::Number, S::Segment)
    # translate and rotate segment to positive Re axis
    d = S.zb - S.za
    s = sign(d)
    ζ = (z - S.za) / s
    return S.za + s * min(max(real(ζ), 0), abs(d))
end

"""
	reflect(z,S::Segment)

Reflect the value `z` across the extension of segment `S` to a line. (For reflection of a segment through a point, use translation and negation.)
"""
reflect(z::Number, S::Segment) = reflect(z, Line(S.za, S.zb))

# Display methods
# COV_EXCL_START
function show(io::IO, S::Segment{T}) where {T}
    print(IOContext(io, :compact => true), "Segment(", point(S, 0), ",", point(S, 1), ")")
end

function show(io::IO, ::MIME"text/plain", S::Segment{T}) where {T}
    print(io, "Segment{$T} in the complex plane:\n   from (", point(S, 0), ") to (", point(S, 1), ")")
end
# COV_EXCL_STOP

plotdata(S::Segment) = [S.za, S.zb]
