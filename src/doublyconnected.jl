# TODO: check that the inner boundary is nested inside the outer
struct InteriorDoublyConnectedRegion{T,S<:Jordan{T},R<:Jordan{T}} <: AbstractDoublyConnectedRegion{T}
    outer::S
    inner::R
end
innerboundary(R::InteriorDoublyConnectedRegion) = [R.inner]    # must be a vector
outerboundary(R::InteriorDoublyConnectedRegion) = R.outer
Base.in(z::Number, R::InteriorDoublyConnectedRegion) = isinside(z, R.outer) && isoutside(z, R.inner)
Base.isfinite(::InteriorDoublyConnectedRegion) = true

# COV_EXCL_START
function Base.show(io::IO, R::InteriorDoublyConnectedRegion)
    outer = isa(R.outer, ClosedCurve) ? "closed curve" : "$(R.outer)"
    inner = isa(R.inner, ClosedCurve) ? "closed curve" : "$(R.inner)"
    if inner == outer == "closed curve"
        print(IOContext(io, :compact => true), "Region between two closed curves")
    else
         print(IOContext(io, :compact => true), "Region between ", outer, " and ", inner)
    end
end

function Base.show(io::IO, ::MIME"text/plain", R::InteriorDoublyConnectedRegion)
    print(io, "Region between $(R.outer) and $(R.inner)")
end

# COV_EXCL_STOP

# TODO: check that the inner boundaries do not nest or intersect
struct ExteriorDoublyConnectedRegion{T,S<:Jordan{T},R<:Jordan{T}} <: AbstractDoublyConnectedRegion{T}
    inner::Tuple{S,R}
end
ExteriorDoublyConnectedRegion(c1::Jordan{T}, c2::Jordan{T}) where T = ExteriorDoublyConnectedRegion((c1, c2))
innerboundary(R::ExteriorDoublyConnectedRegion) = [R.inner...]    # must be a vector
outerboundary(R::ExteriorDoublyConnectedRegion) = nothing
Base.in(z::Number, R::ExteriorDoublyConnectedRegion) = all(isoutside(z, r) for r in R.inner)
Base.isfinite(::ExteriorDoublyConnectedRegion) = false

"""
    between(outer, inner)
Construct the region interior to the closed curve or path `outer` and interior to `inner`.
"""
function between(outer::Jordan{T}, inner::Jordan{T}) where T
    if !all(isinside(outer), discretize(inner, 50)[2])
        if all(isinside(inner), discretize(outer, 50)[2])
            return between(inner, outer)
        end
        throw(ArgumentError("One curve must be inside the other"))
    end
    if isfinite(outer) && isinside(Inf, outer)
        outer = reverse(outer)
    end
    if isfinite(inner) && isoutside(Inf, inner)
        inner = reverse(inner)
    end
    return InteriorDoublyConnectedRegion(outer, inner)
end

#
# Annulus
#

"""
    (type) Annulus

Region between two concentric circles.
"""
struct Annulus{T} <: AbstractDoublyConnectedRegion{T}
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
    return Annulus(Circle(center, outerrad, true), Circle(center, innerrad, false))
end

"""
    modulus(A::Annulus)

Return the conformal modulus of the annulus `A`, defined as the ratio of the inner radius to the outer radius.
"""
modulus(A::Annulus) = A.inner.radius / A.outer.radius
innerboundary(A::Annulus) = [A.inner]    # must be a vector
outerboundary(A::Annulus) = A.outer
Base.in(z::Number, A::Annulus) = isinside(z, A.outer) && isoutside(z, A.inner)
Base.isfinite(::Annulus) = true

# COV_EXCL_START
function Base.show(io::IO, ::MIME"text/plain", R::Annulus)
    print(io, "Annulus centered at ", R.outer.center, " with radii ", R.outer.radius, ", ", R.inner.radius)
end
# COV_EXCL_STOP
