abstract type AbstractParameterizedMap{T} <: Function end

# This is an explicit statement of what is default behavior for a Function. It makes, e.g.,
# point.(c, t) work for parameter arrays, but it disables broadcasting over the
# sides of a Path object. Callers should broadcast over sides(p) explicitly.
Base.broadcastable(::AbstractParameterizedMap) = Ref(x)

# COV_EXCL_START
# Required methods
"""
	point(C::AbstractParameterizedMap, t::Real)
Find the point on curve `C` at parameter value `t`, which should lie in the interval [0,1].
"""
point(c::AbstractParameterizedMap, t::Real) = throw(MethodError(point, (c, t)))

"""
	length(C::AbstractParameterizedMap)
Return the length of the parameter interval for map `C`.
"""
Base.length(::AbstractParameterizedMap) = throw(MethodError(Base.length, (C,)))

"""
	tangent(C::AbstractParameterizedMap, t::Real)
Find the complex number representing the tangent to curve `C` at parameter value `t` in [0,1].
"""
tangent(C::AbstractParameterizedMap, t::Real) = throw(MethodError(tangent, (C, t)))

Base.reverse(C::AbstractParameterizedMap) = throw(MethodError(Base.reverse, (C,)))

arclength(C::AbstractParameterizedMap, args...) = throw(MethodError(arclength, (C, args...)))

"""
	isfinite(C::AbstractParameterizedMap)
Return `true` if the curve is bounded in the complex plane (i.e., does not pass through infinity).
"""
Base.isfinite(C::AbstractParameterizedMap) = throw(MethodError(Base.isfinite, (C,)))
Base.isreal(::AbstractParameterizedMap) = false  # unless detected otherwise

Base.conj(C::AbstractParameterizedMap) = throw(MethodError(Base.conj, (C,)))
Base.:+(C::AbstractParameterizedMap, z::Number) = throw(MethodError(Base.:+, (C, z)))
Base.:-(C::AbstractParameterizedMap) = throw(MethodError(Base.:-, (C,)))
Base.:*(C::AbstractParameterizedMap, z::Number) = throw(MethodError(Base.:*, (C, z)))
Base.inv(C::AbstractParameterizedMap) = throw(MethodError(Base.inv, (C,)))
# COV_EXCL_STOP

# Concrete types must define C + z, C * z, and inv(C) only.
Base.:+(z::Number, C::AbstractParameterizedMap) = C + z
Base.:-(C::AbstractParameterizedMap, z::Number) = C + (-z)
Base.:-(z::Number, C::AbstractParameterizedMap) = z + (-C)
Base.:*(z::Number, C::AbstractParameterizedMap) = C * z
Base.:/(C::AbstractParameterizedMap{T}, z::Number) where {T} = C * (T(1) / z)
Base.:/(z::Number, C::AbstractParameterizedMap) = z * inv(C)

# make callable by name
(C::AbstractParameterizedMap)(t::Real) = point(C, t)

"""
	points(C::AbstractParameterizedMap, t::AbstractArray)

Vectorize the `point` function.
"""
points(C::AbstractParameterizedMap, t::AbstractArray{T}) where {T<:Real} = [point(C, t) for t in t]

@deprecate point(C::AbstractParameterizedMap, t::AbstractArray) points(C, t)

"""
	unittangent(C::AbstractParameterizedMap, t::Real)
Find the complex number representing the unit tangent at parameter value `t`.
"""
unittangent(C::AbstractParameterizedMap, t::Real) = sign(tangent(C, t))

"""
	normal(C::AbstractParameterizedMap,t::Real)
Find the complex number in the direction of the leftward-pointing normal at parameter value `t`.
"""
unitnormal(C::AbstractParameterizedMap, t::Real) = 1im * unittangent(C, t)

@deprecate normal unitnormal
