const AbstractJordan{T} = Union{AbstractClosedCurve{T},AbstractClosedPath{T}}

# get one point inside a closed path
function get_one_inside(C::AbstractJordan{T}) where {T}
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
    @error "Could not find a point inside the curve"
end

abstract type AbstractRegion{T} end

# COV_EXCL_START
# Required methods
boundary(R::AbstractRegion) = @error "No boundary() method defined for type $(typeof(R))"

"""
    in(z::Number,R::AbstractRegion;tol=<default>)
    z ∈ R   (type "\\in" followed by tab)
True if `z` is in the region `R`.
"""
in(z::Number, R::AbstractRegion; tol=nothing) = @error "No in() method defined for type $(typeof(R))"
in(R::AbstractRegion; kw...) = z -> in(z, R; kw...)

"""
    isfinite(R::AbstractRegion)
Return `true` if the region is bounded in the complex plane.
"""
isfinite(R::AbstractRegion) = @error "No isfinite() method defined for type $(typeof(R))"
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

abstract type AbstractConnectedRegion{N,T} <: AbstractRegion{T} end

# Required methods
innerboundary(R::AbstractConnectedRegion) = @error "No innerboundary() method defined for type $(typeof(R))"
outerboundary(R::AbstractConnectedRegion) = @error "No outerboundary() method defined for type $(typeof(R))"
boundary(R::AbstractConnectedRegion) = outerboundary(R), innerboundary(R)

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
function show(io::IO, ::MIME"text/plain", R::AbstractConnectedRegion)
    no = isnothing(outerboundary(R)) ? "no" : ""
    print(io, "Region in the complex plane with $no outer boundary and $(length(innerboundary(R))) inner boundary components")
end

function show(io::IO, R::AbstractConnectedRegion)
    print(io, "Region in the complex plane")
end
# COV_EXCL_STOP

#
# concrete implementations
#

########################
# ExteriorRegion
########################

struct ExteriorRegion{N,T} <: AbstractConnectedRegion{N,T}
    inner::Vector{<:AbstractJordan{T}}
    function ExteriorRegion{N,T}(inner::AbstractVector) where {N,T}
        @assert N == length(inner) "Incorrect connectivity"
        @assert all(c isa AbstractJordan{T} for c in inner) "Boundary components must be closed curves or paths"
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
    try
        T = promote_type(real_type.(inner)...)
        return ExteriorRegion{length(inner),T}(convert(Vector{AbstractJordan{T}}, inner))
    catch
        @error "Could not promote the types of the inner boundaries"
    end
end

isfinite(::ExteriorRegion) = false
in(z::Number, R::ExteriorRegion) = all(isoutside(z, c) for c in R.inner)
innerboundary(R::ExteriorRegion) = R.inner
outerboundary(::ExteriorRegion) = nothing
Base.:+(R::ExteriorRegion, z::Number) = ExteriorRegion(innerboundary(R) .+ z)
Base.:*(R::ExteriorRegion, z::Number) = ExteriorRegion(innerboundary(R) .* z)
Base.:-(R::ExteriorRegion, z::Number) = ExteriorRegion(innerboundary(R) .- z)
Base.:-(R::ExteriorRegion) = ExteriorRegion(-innerboundary(R))

convert_real_type(T::Type{<:AbstractFloat}, R::ExteriorRegion{N,S}) where {N,S} = ExteriorRegion{N,T}(R.inner)
Base.promote_rule(::Type{<:ExteriorRegion{N,T}}, ::Type{<:ExteriorRegion{N,S}}) where {N,T,S} = ExteriorRegion{promote_type(T,S)}

#
# ConnectedRegion
#

"""
    (type) ConnectedRegion{N}
Representation of a `N`-connected region in the extended complex plane.
"""
struct ConnectedRegion{N,T} <: AbstractConnectedRegion{N,T}
    outer::Union{Nothing,AbstractJordan{T}}
    inner::Vector{AbstractJordan{T}}
    function ConnectedRegion{N,T}(outer, inner) where {N,T}
        n = length(inner) + !isnothing(outer)
        @assert N == n "Incorrect connectivity"
        if !isnothing(outer)
            # correct orientation of outer component?
            isin = [isinside(point(c, 0), outer) for c in inner]
            if all(.!isin)
                outer = reverse(outer)
            else
                @assert all(isin) "Inner components appear to be crossing the outer boundary"
            end
        end
        # correct orientations of inner components?
        @assert all(isfinite.(inner)) "Inner boundaries must be finite"
        for c in inner
            if !isfinite(interior(c))
                c = reverse(c)
            end
        end
        new(outer, inner)
    end
end

"""
    ConnectedRegion(outer, inner)
Construct an open connected region by specifying its boundary components. The `outer` boundary could be `nothing` or a closed curve or path. The `inner` boundary should be a vector of one or more nonintersecting closed curves or paths. The defined region is interior to the outer boundary and exterior to all the components of the inner boundary, regardless of the orientations of the given curves.
"""
function ConnectedRegion(
                        outer::Union{Nothing,AbstractJordan{T}},
                        inner::AbstractVector{<:AbstractJordan{T}}
                        ) where T
    n = length(inner)
    if isnothing(outer)
        ExteriorRegion{n,T}(inner)
    else
        ConnectedRegion{n+1,T}(outer, inner)
    end
end

function in(z::Number, R::ConnectedRegion)
    all(isoutside(z, c) for c in R.inner) && isinside(z, R.outer)
end

outerboundary(R::ConnectedRegion) = R.outer
innerboundary(R::ConnectedRegion) = R.inner
#innerboundary(R::ConnectedRegion{2}) = R.inner[1]
convert_real_type(T::Type{<:AbstractFloat}, R::ConnectedRegion{N,S}) where {N,S} = ConnectedRegion{N,T}(R,outer, R.inner)
Base.promote_rule(::Type{ConnectedRegion{N,T}}, ::Type{ConnectedRegion{N,S}}) where {N,T,S} = ConnectedRegion{N,promote_type(T,S)}

#
# special cases
#

"""
    between(outer,inner)
Construct the region interior to the closed curve or path `outer` and interior to `inner`.
"""
function between(outer::AbstractJordan{T}, inner::AbstractJordan{T}) where T
    if isfinite(outer) && isinside(Inf, outer)
        outer = reverse(outer)
    end
    if isfinite(inner) && isoutside(Inf, inner)
        inner = reverse(inner)
    end
    ConnectedRegion{2,T}(outer, [inner])
end


# Annulus

"""
    (type) Annulus
Representation of the region between two circles.
"""
struct Annulus{T} <: AbstractConnectedRegion{2,T}
    outer::Circle{T}
    inner::Circle{T}
    function Annulus{T}(outer::Circle{T}, inner::Circle{T}) where T
        @assert(outer.center ≈ inner.center)
        if isinside(Inf, outer)
            outer = reverse(outer)
        end
        if isoutside(Inf, inner)
            inner = reverse(inner)
        end
        new(outer, inner)
    end
end

"""
    Annulus(radouter, radinner)
    Annulus(radouter, radinner, center)
Construct a concentric annulus of outer radius `radouter` and inner radius `radinner` centered at `center`. If the center is not given, the origin is used.
"""
function Annulus(outer::Circle{T}, inner::Circle{S}) where {T,S}
    R = promote_type(T,S)
    return Annulus{R}(convert(Circle{R}, outer), convert(Circle{R}, inner))
end

function Annulus(outerrad::Real, innerrad::Real, center::Number=0)
    @assert outerrad > innerrad > 0
    Annulus(Circle(center, outerrad, true), Circle(center, innerrad, false))
end

modulus(A::Annulus) = A.inner.radius / A.outer.radius
innerboundary(A::Annulus) = [A.inner]    # must be a vector
outerboundary(A::Annulus) = A.outer
in(z::Number, A::Annulus) = isinside(z, A.outer) && isoutside(z, A.inner)
isfinite(::Annulus) = true

# COV_EXCL_START
function show(io::IO, ::MIME"text/plain", R::Annulus)
    print(io, "Annulus in the complex plane:\n")
    print(io, "   centered at ", R.outer.center, " with distances from ", R.inner.radius, " to ", R.outer.radius)
end
# COV_EXCL_STOP
