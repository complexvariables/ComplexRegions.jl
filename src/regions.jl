const Jordan{T} = Union{AbstractClosedCurve{T},AbstractClosedPath{T}}

# get one point inside a closed path
function get_one_inside(C::Jordan{T}) where {T}
    zc = mean(discretize(C, ds=1/T(100)))
    if isinside(zc, C)   # ignores the orientation
        return zc
    elseif isinside(0, C)  # allows manual control
        return 0
    else
        # Try midpoints between pairs of boundary points
        t = range(T(0), T(length(C)), 11)
        for a in t, b in t .+ T(1)/20
            zc = mean(C.([a, b]))
            if isinside(zc, C)
                return zc
            end
        end
    end
    error("Could not find a point inside the curve")
end

abstract type AbstractRegion{T} end

# COV_EXCL_START
# Required methods
boundary(R::AbstractRegion) = throw(MethodError(boundary, (R,)))

"""
    in(z::Number,R::AbstractRegion;tol=<default>)
    z ∈ R   (type "\\in" followed by tab)
True if `z` is in the region `R`.
"""
in(z::Number, R::AbstractRegion; tol=nothing) = throw(MethodError(in, (z, R)))
in(R::AbstractRegion; kw...) = z -> in(z, R; kw...)

"""
    isfinite(R::AbstractRegion)
Return `true` if the region is bounded in the complex plane.
"""
isfinite(R::AbstractRegion) = throw(MethodError(isfinite, (R,)))
# COV_EXCL_STOP

# Default implementations
"""
    (type) RegionIntersection
Representation of the intersection of two regions.
"""
struct RegionIntersection{T} <: AbstractRegion{T}
    one::AbstractRegion{T}
    two::AbstractRegion{T}
    function RegionIntersection(one::AbstractRegion{R}, two::AbstractRegion{S}) where {R,S}
        # one, two = promote(one, two)
        new{promote_type(R,S)}(one, two)
    end
end
in(z::Number, R::RegionIntersection) = in(z, R.one) && in(z, R.two)

"""
    (type) RegionUnion
Representation of the union of two regions.
"""
struct RegionUnion{T} <: AbstractRegion{T}
    one::AbstractRegion{T}
    two::AbstractRegion{T}
end
in(z::Number, R::RegionUnion) = in(z, R.one) || in(z, R.two)

"""
    intersect(R1::AbstractRegion,R2::AbstractRegion)
    R1 ∩ R2    (type "\\cap" followed by tab key)
Create the region that is the intersection of `R1` and `R2`.
"""
function intersect(R1::AbstractRegion{T}, R2::AbstractRegion{S}) where {S,T}
    return RegionIntersection(R1, R2)
end

"""
    union(R1::AbstractRegion,R2::AbstractRegion)
    R1 ∪ R2    (type "\\cup" followed by tab key)
Create the region that is the union of `R1` and `R2`.
"""
function union(R1::AbstractRegion{T}, R2::AbstractRegion{S}) where {S,T}
    return RegionUnion{promote_type(S,T)}(R1, R2)
end

#############################
# AbstractConnectedRegion
#############################

abstract type AbstractConnectedRegion{T} <: AbstractRegion{T} end
abstract type AbstractSimplyConnectedRegion{T} <: AbstractConnectedRegion{T} end
abstract type AbstractDoublyConnectedRegion{T} <: AbstractConnectedRegion{T} end

# Required methods
innerboundary(R::AbstractConnectedRegion) = throw(MethodError(innerboundary, (R,)))
outerboundary(R::AbstractConnectedRegion) = throw(MethodError(outerboundary, (R,)))
boundary(R::AbstractConnectedRegion) = outerboundary(R), innerboundary(R)

"""
    connectivity(R::AbstractConnectedRegion)
Return the connectivity of `R` (the number of connected boundary components).
"""
connectivity(R::AbstractSimplyConnectedRegion) = 1
connectivity(R::AbstractDoublyConnectedRegion) = 2
function connectivity(R::AbstractConnectedRegion)
    n = length(innerboundary(R))
    isnothing(outerboundary(R)) ? n : n + 1
end

# Default implementations

Base.:+(R::AbstractConnectedRegion, z::Number) = typeof(R)(outerboundary(R) + z, innerboundary(R) .+ z)
Base.:+(z::Number, R::AbstractConnectedRegion) = +(R, z)

Base.:-(R::AbstractConnectedRegion) = typeof(R)(-outerboundary(R), map(-, innerboundary(R)))
Base.:-(R::AbstractConnectedRegion, z::Number) = +(R, -z)
Base.:-(z::Number, R::AbstractConnectedRegion) = +(z, -R)

Base.:*(R::AbstractConnectedRegion, z::Number) = typeof(R)(outerboundary(R) * z, innerboundary(R) .* z)
Base.:*(z::Number, R::AbstractConnectedRegion) = R * z

Base.:/(R::AbstractConnectedRegion, z::Number) = *(R, 1 / z)
#/(z::Number,R::AbstractConnectedRegion) = z*inv(R)
#inv(p::AbstractConnectedRegion) = typeof(R)([inv(c) for c in curves(p)])

# COV_EXCL_START
function Base.show(io::IO, ::MIME"text/plain", R::AbstractConnectedRegion)
    no = isnothing(outerboundary(R)) ? "no " : ""
    n = length(innerboundary(R))
    print(io, "Region in the complex plane with $(no)outer boundary and $(n) inner boundary component")
    if n > 1
        print(io, "s")
    end
end

function Base.show(io::IO, R::AbstractConnectedRegion)
    print(io, "Region in the complex plane")
end
# COV_EXCL_STOP

########################
# generic ExteriorRegion
########################

struct ExteriorRegion{T} <: AbstractConnectedRegion{T}
    inner::Vector{<:Jordan{T}}
    function ExteriorRegion{T}(inner::AbstractVector) where {T}
        @assert all(c isa Jordan{T} for c in inner) "Boundary components must be closed curves or paths"
        @assert all(isfinite.(inner)) "Inner boundaries must be finite"
        # Correct the orientations of inner components (region is on the left)
        b = ClosedPath.(copy(inner))
        for k in eachindex(b)
            if winding(b[k], get_one_inside(b[k])) > 0
                b[k] = reverse(b[k])
            end
        end
        new(b)
    end
end

function ExteriorRegion(inner::AbstractVector)
    T = promote_type(real_type.(inner)...)
    return ExteriorRegion{T}(convert(Vector{Jordan{T}}, inner))
end

Base.isfinite(::ExteriorRegion) = false
Base.in(z::Number, R::ExteriorRegion) = all(isoutside(z, c) for c in R.inner)
innerboundary(R::ExteriorRegion) = R.inner
outerboundary(::ExteriorRegion) = nothing
Base.:+(R::ExteriorRegion, z::Number) = ExteriorRegion(innerboundary(R) .+ z)
Base.:*(R::ExteriorRegion, z::Number) = ExteriorRegion(innerboundary(R) .* z)
Base.:-(R::ExteriorRegion, z::Number) = ExteriorRegion(innerboundary(R) .- z)
Base.:-(R::ExteriorRegion) = ExteriorRegion(-innerboundary(R))

convert_real_type(T::Type{<:AbstractFloat}, R::ExteriorRegion{S}) where {S} = ExteriorRegion{T}(R.inner)
Base.promote_rule(::Type{<:ExteriorRegion{T}}, ::Type{<:ExteriorRegion{S}}) where {T,S} = ExteriorRegion{promote_type(T,S)}

########################
# generic InteriorRegion
########################

"""
    (type) InteriorRegion{N}
Representation of a `N`-connected region in the extended complex plane.
"""
struct InteriorRegion{T} <: AbstractConnectedRegion{T}
    outer::Jordan{T}
    inner::Vector{<:Jordan{T}}
    function InteriorRegion(outer::Jordan, inner::AbstractVector{<:Jordan})
        T = promote_type(real_type(outer), real_type.(inner)...)
        return InteriorRegion{T}(outer, inner)
    end
    function InteriorRegion{T}(outer, inner) where {T<:AbstractFloat}
        # correct orientation of outer component?
        isin = [isinside(point(c, 0), outer) for c in inner]
        if all(.!isin)
            outer = reverse(outer)
        else
            @assert all(isin) "Inner components appear to be crossing the outer boundary"
        end
        # correct orientations of inner components?
        @assert all(isfinite.(inner)) "Inner boundaries must be finite"
        for c in inner
            if !isfinite(interior(c))
                c = reverse(c)
            end
        end
        new{T}(outer, inner)
    end
end

"""
    connected_region(outer, inner)
Construct an open connected region by specifying its boundary components. The `outer` boundary could be `nothing` or a closed curve or path. The `inner` boundary should be a vector of one or more nonintersecting closed curves or paths. The defined region is interior to the outer boundary and exterior to all the components of the inner boundary, regardless of the orientations of the given curves.
"""
function connected_region(inner::AbstractVector{<:Jordan{T}}) where T
    ExteriorRegion{T}(inner)
end

function connected_region(
                        outer::Jordan{T},
                        inner::AbstractVector{<:Jordan{T}}
                        ) where T
    InteriorRegion{T}(outer, inner)
end


function Base.in(z::Number, R::InteriorRegion)
    all(isoutside(z, c) for c in R.inner) && isinside(z, R.outer)
end

outerboundary(R::InteriorRegion) = R.outer
innerboundary(R::InteriorRegion) = R.inner
convert_real_type(T::Type{<:AbstractFloat}, R::InteriorRegion{S}) where {S} = InteriorRegion{T}(R.outer, R.inner)
Base.promote_rule(::Type{InteriorRegion{T}}, ::Type{InteriorRegion{S}}) where {T,S} = InteriorRegion{promote_type(T,S)}

@deprecate InteriorConnectedRegion InteriorRegion