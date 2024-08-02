abstract type AbstractPath{T} end
const PathLike{T} = Union{AbstractPath{T},AbstractVector{<:AbstractCurve{T}}}

# Required methods
"""
	curves(P::AbstractPath)
Return an array of the curves that make up the path `P`.
"""
curves(p::AbstractPath)::AbstractVector = @error "No curves() method defined for type $(typeof(p))"

# Default implementations
"""
	curve(P::AbstractPath,k::Integer)
Return the `k`th curve in the path `P`.
"""
curve(p::AbstractPath, k::Integer) = curves(p)[k]

"""
	vertex(P::AbstractPath,k::Integer)
Return the `k`th vertex of the path `P`.
"""
function vertex(P::AbstractPath, k::Integer)
    C = curves(P)
    n = length(C)
    if 1 ≤ k ≤ n
        point(C[k], 0)
    elseif k == n + 1
        point(C[n], 1)
    else
        throw(BoundsError(P, k))
    end
end

"""
	vertices(P::AbstractPath)
Return an array of the vertices (endpoints of the curves) of the path `P`. The length is one greater than the number of curves in `P`.
"""
function vertices(P::AbstractPath)
    [vertex(P, k) for k = 1:length(P)+1]
end

"""
	isfinite(P::AbstractPath)
Return `true` if the path is bounded in the complex plane (i.e., does not pass through infinity).
"""
isfinite(p::AbstractPath) = all(isfinite(s) for s in curves(p))

"""
	isreal(P::AbstractPath)
Return `true` if the path is entirely on the real axis.
"""
Base.isreal(p::AbstractPath) = all(isreal(s) for s in curves(p))

"""
	real_type(::AbstractPath)
Return the type of the real part of the curve's point function.
"""
real_type(::AbstractPath{T}) where T = T

# iteration interface
eltype(::Type{AbstractPath{T}}) where T = AbstractCurve{T}
length(p::AbstractPath) = length(curves(p))
getindex(p::AbstractPath, k) = curve(p, k)
iterate(p::AbstractPath, state=1) = state > length(curves(p)) ? nothing : (p[state], state + 1)

"""
	point(P::AbstractPath, t::Real)
	P(t)
Compute the point along path `P` at parameter value `t`. Values of `t` in [k,k+1] correspond to values in [0,1] along curve k of the path, for k = 1,2,...,length(P)-1.

	point(P::AbstractPath, t::AbstractVector)
Vectorize the `point` method for path `P`.
"""
point(p::AbstractPath, t::AbstractArray{<:Real}) = [point(p, t) for t in t]
(p::AbstractPath)(t) = point(p, t)

# This determines how to parse a parameter of a path. Overloaded later for the case of a closed path.
function sideargs(p::AbstractPath, t)
    n = length(p)
    if (t < 0) || (t > n)
        throw(BoundsError(p, t))
    end
    if t == n
        return n, 1
    else
        return 1 + floor(Int, t), t % 1
    end
end

function point(p::AbstractPath{T}, t::Real) where T
    k, s = sideargs(p, T(t))
    point(curve(p, k), s)
end

"""
	tangent(P::AbstractPath,t::Real)
Compute the complex-valued tangent along path `P` at parameter value `t`. Values of `t` in [k,k+1] correspond to values in [0,1] along curve k of the path, for k = 1,2,...,length(P)-1. The result is not well-defined at an integer value of `t`.
"""
function tangent(p::AbstractPath{T}, t::Real) where T
    k, s = sideargs(p, T(t))
    tangent(curve(p, k), s)
end

"""
	unittangent(P::AbstractPath,t::Real)
Compute the complex-valued unit tangent along path `P` at parameter value `t`. Values of `t` in [k,k+1] correspond to values in [0,1] along curve k of the path, for k = 1,2,...,length(P)-1. The result is not well-defined at an integer value of `t`.
"""
function unittangent(p::AbstractPath{T}, t::Real) where T
    k, s = sideargs(p, T(t))
    unittangent(curve(p, k), s)
end

"""
	normal(P::AbstractPath,t::Real)
Compute a complex-valued normal to path `P` at parameter value `t`. Values of `t` in [k,k+1] correspond to values in [0,1] along curve k of the path, for k = 1,2,...,length(P)-1. The result is not well-defined at an integer value of `t`.
"""
function normal(p::AbstractPath{T}, t::Real) where T
    k, s = sideargs(p, T(t))
    normal(curve(p, k), s)
end

"""
	angles(P::AbstractPath)
Return a vecrtor of the interior angles at the vertices of the path `P`. The length is one greater than the number of curves in `P`, and the first and last values are `NaN`.
"""
function angles(P::AbstractPath{T}) where T
    τ = nothing
    try
        τ = tangent(P, 0)
    catch
        @error "Path must have a tangent defined"
    end
    m = length(P)
    θ = Vector{real_type(τ)}(undef, m + 1)
    θ[1] = θ[m+1] = NaN
    c = curves(P)
    for n in 1:m-1
        τminus = -tangent(c[n], 1)
        τplus = tangent(c[n+1], 0)
        θ[n+1] = mod2pi(angle(τminus / τplus))
    end
    return θ
end

conj(p::AbstractPath) = typeof(p)(conj.(curves(p)))

reverse(p::AbstractPath) = typeof(p)(reverse(reverse.(curves(p))))
isclosed(p::AbstractPath) = isa(p, AbstractClosedPath)

function Base.:+(p::AbstractPath{T}, z::Number) where T
	return typeof(p)([c + convert_real_type(T, z) for c in curves(p)])
end
Base.:+(z::Number, p::AbstractPath) = p + z

Base.:-(p::AbstractPath) = typeof(p)([-c for c in curves(p)])
Base.:-(p::AbstractPath, z::Number) = p + (-z)
Base.:-(z::Number, p::AbstractPath) = (-p) + z

Base.:*(p::AbstractPath, z::Number) = typeof(p)([c * z for c in curves(p)])
Base.:*(z::Number, p::AbstractPath) = typeof(p)([z * c for c in curves(p)])

Base.:/(p::AbstractPath, z::Number) = typeof(p)([c / z for c in curves(p)])
Base.:/(z::Number, p::AbstractPath) = z * inv(p)
inv(p::AbstractPath) = typeof(p)([inv(c) for c in curves(p)])

"""
	isapprox(P1::AbstractPath,R2::AbstractPath)
	P1 ≈ P2       (type "\\approx" followed by tab key)
Determine whether `P1` and `P2` represent the same path, up to tolerance `tol`, irrespective of the parameterization of its curves.
"""
function isapprox(P1::AbstractPath, P2::AbstractPath)
    if length(P1) != length(P2)
        return false
    else
        c1, c2 = curves(P1), curves(P2)
        all(isapprox(c1[k], c2[k]) for k in eachindex(c1))
    end
end
isapprox(::AbstractCurve, ::AbstractPath; kw...) = false
isapprox(::AbstractPath, ::AbstractCurve; kw...) = false

"""
	dist(z,P::AbstractPath)
Find the distance from the path `P` to the point `z`.
"""
dist(z::Number, P::AbstractPath) = minimum(dist(z, C) for C in P)

"""
	closest(z,P::AbstractPath)
Find the point on the path `P` that lies closest to `z`.
"""
function closest(z::Number, P::AbstractPath)
    k = argmin([dist(z, s) for s in curves(P)])
    closest(z, side(P, k))
end

intersect(P::AbstractPath, C::AbstractCurve) = ∪([intersect(s, C) for s in curves(P)]...)
intersect(C::AbstractCurve, P::AbstractPath) = intersect(P, C)
intersect(P1::AbstractPath, P2::AbstractPath) = ∪([intersect(P1, s) for s in curves(P2)]...)

function show(io::IO, P::AbstractPath)
    str = length(P) == 1 ? " curve" : " curves"
    print(IOContext(io, :compact => true), typeof(P), " with ", length(P), str)
end
function show(io::IO, ::MIME"text/plain", P::AbstractPath)
    str = length(P) == 1 ? " curve" : " curves"
    print(io, typeof(P), " with ", length(P), str)
end

###############
# ClosedPath
###############

abstract type AbstractClosedPath{T} <: AbstractPath{T} end

"""
	curve(P::AbstractClosedPath,k::Integer)
Return the `k`th curve in the path `P`. The index is applied circularly; e.g, if the closed path has n curves, then ...,1-n,1,1+n,... all refer to the first curve.
"""
curve(p::AbstractClosedPath, k::Integer) = curves(p)[mod(k - 1, length(p))+1]

"""
	vertex(P::AbstractPath,k::Integer)
Return the `k`th vertex of the path `P`. The index is applied circularly; e.g, if the closed path has n curves, then ...,1-n,1,1+n,... all refer to the first vertex.
"""
vertex(P::AbstractClosedPath, k::Integer) = point(curve(P, k), 0)

"""
	vertices(P::AbstractClosedPath)
Return an array of the unique vertices (endpoints of the curves) of the closed path `P`. The length is equal the number of curves in `P`, i.e., the first/last vertex is not duplicated.
"""
vertices(P::AbstractClosedPath) = [vertex(P, k) for k = 1:length(P)]

"""
	angles(P::AbstractPath)
Return a vecrtor of the interior angles at the vertices of the path `P`. The length is one greater than the number of curves in `P`, and the first and last values are `NaN`.
"""
function angles(P::AbstractClosedPath{T}) where T
    τ = nothing
    try
        τ = tangent(P, 0)
    catch
        @error "Path must have a tangent defined"
    end
    m = length(P)
    θ = similar([real(τ)], m)
    c = CircularVector(curves(P))
    for n in 0:m-1
        τminus = -tangent(c[n], 1)
        τplus = tangent(c[n+1], 0)
        θ[n+1] = mod2pi(angle(τminus / τplus))
    end
    return θ
end

function sideargs(p::AbstractClosedPath{T}, t) where T
    n = length(p)
    return 1 + mod(floor(Int, T(t)), n), mod(T(t), 1)
end

function winding(P::AbstractClosedPath, z::Number)
    # Integrate around the boundary
    w = 0
    for s in P
        f = t -> imag(tangent(s, t) / (point(s, t) - z))
        w += intadapt(f, 0, 1, 1e-3)
    end
    return round(Int, w / (2π))
end

isinside(z::Number, P::AbstractClosedPath) = winding(P, z) != 0
isoutside(z::Number, P::AbstractClosedPath) = winding(P, z) == 0

#
# Concrete implementations
#

# Path
"""
	(type) Path
Generic implementation of an `AbstractPath`.
"""
struct Path{T} <: AbstractPath{T}
    curve::Vector{AbstractCurve{T}}
    function Path{T}(c::AbstractVector{<:AbstractCurve{T}}; tol::Real=tolerance(T)) where {T}
        n = length(c)
        for k = 1:n-1
            @assert isapprox(point(c[k], 1), point(c[k+1], 0), rtol=tol, atol=tol) "Curve endpoints do not match for pieces $(k) and $(k+1)"
        end
        new{T}(c)
    end
end
"""
	Path(c::AbstractVector; tol=<default>)
Given a vector `c` of curves, construct a path. The path is checked for continuity (to tolerance `tol`) at the interior vertices.
"""
Path(c::AbstractVector{<:AbstractCurve{T}}) where {T} = Path{T}(c)
function Path(p::AbstractVector)
    try
        T = promote_type(real_type.(p)...)
        return Path{T}(convert(Vector{Curve{T}}, p))
    catch
        @error "Vector must contain Curves"
    end
end
Path(c::AbstractCurve{T}) where {T} = Path{T}([c])

curves(p::Path) = p.curve
"""
	arclength(P::Path)
Compute the arclength of the path `P`.
"""
arclength(p::Path) = sum(arclength(c) for c in p)

"""
	convert_real_type(T::Type{<:AbstractFloat}, Union{Path{S},ClosedPath{S}})
Convert the floating-point type of the point and tangent functions of a `Path` or `ClosedPath` object.
"""
function convert_real_type(T::Type{<:Real}, P::Path{S}) where S
    return Path{T}(convert_real_type.(Ref(T), curves(P)))
end

# ClosedPath
"""
	(type) ClosedPath
Generic implementation of an `AbstractClosedPath`.
"""
struct ClosedPath{T} <: AbstractClosedPath{T}
    curve::Vector{AbstractCurve{T}}
    function ClosedPath{T}(p::AbstractVector{<:AbstractCurve{T}}; tol=tolerance(T)) where {T}
        q = Path{T}(p)
        zi, zf = point(q, 0), point(q, length(q))
        @assert isapprox(zi, zf, rtol=tol, atol=tol) || (isinf(zi) && isinf(zf)) "Path endpoints do not match"
        new{T}(p)
    end
end
"""
	ClosedPath(c::AbstractVector; tol=<default>)
	ClosedPath(P::Path; tol=<default>)
Given a vector `c` of curves, or an existing path, construct a closed path. The path is checked for continuity (to tolerance `tol`) at all of the vertices.
"""
ClosedPath(p::AbstractVector{<:AbstractCurve{T}}) where {T} = ClosedPath{T}(p)
function ClosedPath(p::AbstractVector)
    try
        T = promote_type(real_type.(p)...)
        return ClosedPath{T}(convert(Vector{Curve{T}}, p))
    catch
        @error "Vector must contain Curves"
    end
end
ClosedPath(c::AbstractCurve{T}) where {T} = ClosedPath{T}([c])
ClosedPath(c::AbstractClosedPath{T}) where {T} = c
ClosedPath(p::Path; kw...) = ClosedPath(p.curve; kw...)

Base.convert(::Type{ClosedPath}, c::AbstractClosedCurve) = ClosedPath([c])

function convert_real_type(T::Type{<:Real}, P::ClosedPath{S}) where S
    return ClosedPath{T}(convert_real_type.(Ref(T), curves(P)))
end

curves(p::ClosedPath) = p.curve
arclength(p::ClosedPath) = sum(arclength(c) for c in p)

# Find a circle that fully encloses all the finite vertices and some points of a path.
function enclosing_circle(p::AbstractPath, expansion=2)
    z = discretize(p, ds=0.02)
    return Circle(enclosing_circle(filter(isfinite, z), expansion)...)
end

function enclosing_circle(p::AbstractVector{<:AbstractPath}, expansion=2)
    z = []
    map(p) do p
        z = vcat(z, discretize(p, ds=0.02))
    end
    return Circle(enclosing_circle(filter(isfinite, z), expansion)...)
end

include("polygons.jl")
