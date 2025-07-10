#
# abstract interfaces
#

abstract type AbstractCurve{T} <: Function end

# Required methods
"""
	point(C::AbstractCurve,t::Real)
Find the point on curve `C` at parameter value `t`, which should lie in the interval [0,1].
"""
point(c::AbstractCurve, ::Real) = @error "No point() method defined for type $(typeof(c))"

"""
	tangent(C::AbstractCurve,t::Real)
Find the complex number representing the tangent to curve `C` at parameter value `t` in [0,1].
"""
tangent(C::AbstractCurve, ::Real) = @error "No tangent() method defined for type $(typeof(C))"

Base.reverse(C::AbstractCurve) = @error "No reverse() method defined for type $(typeof(C))"

"""
	isfinite(C::AbstractCurve)
Return `true` if the curve is bounded in the complex plane (i.e., does not pass through infinity).
"""
Base.isfinite(C::AbstractCurve) = @error "No isfinite() method defined for type $(typeof(C))"
Base.isreal(::AbstractCurve) = false  # unless detected otherwise

Base.conj(C::AbstractCurve) = @error "No conj() method defined for type $(typeof(C))"
Base.:+(C::AbstractCurve, z::Number) = @error "No addition method defined for type $(typeof(C))"
Base.:-(C::AbstractCurve) = @error "No negation method defined for type $(typeof(C))"
Base.:*(C::AbstractCurve, z::Number) = @error "No multiplication method defined for type $(typeof(C))"
Base.inv(C::AbstractCurve) = @error "No inversion method defined for type $(typeof(C))"

# Default implementations
Base.length(::AbstractCurve) = 1  # length of parameter interval

"""
	real_type(::AbstractCurve)
Return the type of the real part of the curve's point function.
"""
real_type(::AbstractCurve{T}) where T = T

"""
	point(C::AbstractCurve, t::AbstractArray)

Vectorize the `point` function for curve `C`.
"""
point(c::AbstractCurve, t::AbstractArray{T}) where {T<:Real} = [point(c, t) for t in t]

# callable by name
(c::AbstractCurve)(t::Real) = point(c, t)

"""
	unittangent(C::AbstractCurve, t::Real)
Find the complex number representing the unit tangent to curve `C` at parameter value `t` in [0,1]. For Lines, Segments, and Rays, the `t` argument is optional.
"""
unittangent(C::AbstractCurve, t::Real) = sign(tangent(C, t))

"""
	normal(C::AbstractCurve,t::Real)
Find the unit complex number in the direction of the leftward-pointing normal to curve `C` at parameter value `t` in [0,1].
"""
normal(c::AbstractCurve, t::Real) = 1im * unittangent(c, t)

# Concrete types must define C + z, C * z, and inv(C) only.
Base.:+(z::Number, C::AbstractCurve) = C + z
Base.:-(C::AbstractCurve, z::Number) = C + (-z)
Base.:-(z::Number, C::AbstractCurve) = z + (-C)
Base.:*(z::Number, C::AbstractCurve) = C * z
Base.:/(C::AbstractCurve, z::Number) = C * (1 / z)
Base.:/(z::Number, C::AbstractCurve) = z * inv(C)

function arclength(C::AbstractCurve{T}, part=[0, 1]) where T
    f = t -> abs(tangent(C, T(t)))
    intadapt(f, part..., tolerance(T))
end

isclosed(c::AbstractCurve) = isa(c, AbstractClosedCurve)

# COV_EXCL_START
Base.show(io::IO, C::AbstractCurve) = print(io, "Complex-valued $(typeof(C))")
Base.show(io::IO, ::MIME"text/plain", C::AbstractCurve) = print(io, "Complex-valued $(typeof(C))")
# COV_EXCL_END

# By default, curves are not equal to one another, but curves of the same type can override this.
isapprox(::AbstractCurve, ::AbstractCurve; kw...) = false

# AbstractClosedCurve
abstract type AbstractClosedCurve{T} <: AbstractCurve{T} end

# Default implementations
function winding(C::AbstractClosedCurve{T}, z::Number) where T
	Tz = convert_real_type(T, z)
    # Integrate around the curve
    f = t -> imag(tangent(C, T(t)) / (point(C, T(t)) - Tz))
    try
		w = intadapt(f, 0, 1, 1e-4)
    	return round(Int, w / (2π))
	catch
		@error "Unable to determine winding number at $(z)"
	end
end
winding(C::AbstractClosedCurve) = z -> winding(C, z)

isinside(z::Number, C::AbstractClosedCurve) = winding(C, z) != 0
isinside(C::AbstractClosedCurve) = z -> isinside(z, C)
isoutside(z::Number, C::AbstractClosedCurve) = winding(C, z) == 0
isoutside(C::AbstractClosedCurve) = z -> isoutside(z, C)

##################
# generic Curve
##################

"""
(type) Smooth curve defined by an explicit function of a real paramerter in [0,1].
"""
struct Curve{T} <: AbstractCurve{T}
    point::Function
    tangent::Function
end

"""
	Curve(f)
	Curve(f, a, b)
Construct a `Curve` object from the complex-valued function `f` accepting an argument in the interval [0,1]. If `a` and `b` are given, they are the limits of the parameter in the call to the supplied `f`. However, the resulting object will be defined on [0,1], which is internally scaled to [a,b].

	Curve(f, df[, a, b])
Construct a curve with point location and tangent given by the complex-valued functions `f` and `df`, respectively, optionally with given limits on the parameter.
"""

Curve{T}(f::Function, a::Real, b::Real) where T = Curve{T}(f, t -> ForwardDiff.derivative(f, T(t)), a, b)
function Curve{T}(f::Function, df::Function, a::Real, b::Real) where T<:AbstractFloat
    return Curve{T}(t -> f(scaleto(T, a, b, t)), t -> df(scaleto(T, a, b, t)) * (b - a))
end

# Constructors that try to determine underlying Real type automatically
function Curve(
	f::Function,
	df::Union{Function,Nothing}=nothing,
	)
	T = real_type(float(f(0.1234)))
	# didn't know T at call time, so manually handle default argument for df
	if isnothing(df)
		return Curve{T}(f, t -> ForwardDiff.derivative(f, T(t)))
	else
		return Curve{T}(f, df)
	end
end

function Curve(f::Function, a::Real, b::Real)
	T = promote_type(real_type(float(a)), real_type(float(b)))
	Curve{T}(f, t -> ForwardDiff.derivative(f, T(t)), T(a), T(b))
end

function Curve(f::Function, df::Function, a::Real, b::Real)
	T = promote_type(real_type(float(a)), real_type(float(b)))
	Curve{T}(f, df, T(a), T(b))
end

"""
	convert_real_type(T::Type{<:AbstractFloat}, ::Union{Curve{S},ClosedCurve{S}})
Convert the floating-point type of the point and tangent functions of a `Curve` or `ClosedCurve` object.
"""
convert_real_type(T::Type{<:AbstractFloat}, C::Curve{S}) where S = Curve{T}(C.point, C.tangent)
Base.promote_rule(::Type{<:Curve{T}}, ::Type{<:Curve{S}}) where {T,S} = Curve{promote_type(T,S)}

# Required methods
point(C::Curve{T}, t::Real) where T = C.point(T(t))
tangent(C::Curve{T}, t::Real) where T = C.tangent(T(t))
Base.conj(C::Curve{T}) where T = Curve{T}(t -> conj(C.point(t)), t -> conj(C.tangent(t)))
Base.reverse(C::Curve{T}) where T = Curve{T}(t -> C.point(1 - T(t)), t -> -C.tangent(1 - T(t)))
Base.isfinite(C::Curve) = true    # unbounded curves must be explicitly typed

function Base.:+(C::Curve{T}, z::Number) where T
	return Curve{T}(t -> C.point(t) + convert_real_type(T, z), C.tangent)
end

function Base.:*(C::Curve{T}, z::Number) where T
	Tz = convert_real_type(T, z)
	return Curve{T}(t -> C.point(t) * Tz, t -> C.tangent(t) * Tz)
end

Base.:-(C::Curve{T}) where T = Curve{T}(t -> -C.point(t), t -> -C.tangent(t))
Base.inv(C::Curve{T}) where T = Curve{T}(t -> 1 / C.point(t), t -> -C.tangent(t) / C.point(t)^2)

"""
	plotdata(C::AbstractCurve)

Compute points along the curve `C` suitable to make a nice plot of it.
"""
plotdata(C::ComplexRegions.AbstractCurve) = adaptpoints(t -> point(C,t), t -> unittangent(C,t), 0, 1)

#
# generic ClosedCurve
#

"""
(type) Smooth closed curve defined by an explicit function of a real paramerter in [0,1].
"""
struct ClosedCurve{T} <: AbstractClosedCurve{T}
    curve::Curve{T}
    function ClosedCurve{T}(c::Curve; tol=tolerance(T)) where T<:AbstractFloat
        @assert isapprox(point(c, T(0)), point(c, T(1)); rtol=tol, atol=tol) "Curve does not close"
        new(convert_real_type(T, c))
    end
	ClosedCurve{T}(f::Function, args...) where T = ClosedCurve{T}(Curve{T}(f, args...))
end

"""
	ClosedCurve(f; tol=<default>)
	ClosedCurve(f,a,b; tol=<default)
Construct a `ClosedCurve` object from the complex-valued function `point` accepting an argument in the interval [0,1]. The constructor checks whether `f(0)≈f(1)` to tolerance `tol`. If `a` and `b` are given, they are the limits of the parameter in the call to the supplied `f`. However, the resulting object will be defined on [0,1], which is internally scaled to [a,b].

	ClosedCurve(f,df[,a,b]; tol=<default>)
Construct a closed curve with point location and tangent given by the complex-valued functions `f` and `df`, respectively, optionally with given limits on the parameter.
"""
ClosedCurve(f,df=t -> ForwardDiff.derivative(f,t);kw...) = ClosedCurve(Curve(f,df;kw...))
ClosedCurve(f,a::Real,b::Real;kw...) = ClosedCurve(Curve(f,a,b;kw...))
ClosedCurve(f,df,a::Real,b::Real;kw...) = ClosedCurve(Curve(f,df,a,b;kw...))

ClosedCurve(c::Curve{T}, args...) where T = ClosedCurve{T}(c, args...)
function ClosedCurve(f::Function, args...)
	c = Curve(f, args...)
	return ClosedCurve{real_type(c)}(c)
end

point(C::ClosedCurve, t::Real) = point(C.curve, t)
tangent(C::ClosedCurve, t::Real) = tangent(C.curve, t)
Base.conj(C::ClosedCurve) = ClosedCurve(conj(C.curve))
Base.reverse(C::ClosedCurve) = ClosedCurve(reverse(C.curve))
Base.isfinite(C::ClosedCurve) = isfinite(C.curve)
Base.:+(C::ClosedCurve, z::Number) = ClosedCurve(C.curve + z)
Base.:-(C::ClosedCurve) = ClosedCurve(-C.curve)
Base.:*(C::ClosedCurve, z::Number) = ClosedCurve(C.curve * z)
Base.inv(C::ClosedCurve) = ClosedCurve(inv(C.curve))

convert_real_type(T::Type{<:Real}, C::ClosedCurve{S}) where S = ClosedCurve{T}(C.point, C.tangent)
Base.promote_rule(::Type{<:ClosedCurve{T}}, ::Type{<:ClosedCurve{S}}) where {T,S} = ClosedCurve{promote_type(T,S)}

include("lines.jl")
include("rays.jl")
include("segments.jl")
include("circles.jl")
include("arcs.jl")
include("intersections.jl")
