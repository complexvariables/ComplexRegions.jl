# Regions

A **region** is the open set enclosed by a closed curve or closed path.

## [Abstract interface](@id interface_regions)

An `AbstractRegion` is expected to provide the following methods. Here `R` is an `AbstractRegion` and `z` is a number.

`boundary(R)`\
Return the boundary of a region. The number and type of outputs depend on the implementation of the concrete type.

`in(z,R)`\
True if `z` is in the region `R`.

There are default implementations of the following methods.

`isfinite(R)`\
True if the region is bounded in the complex plane.

`intersect(R1,R2)` or `R1 ∩ R2`\
Construct the intersection of two regions, returning a value of subtype `RegionIntersection`. (Currently a stub for further development.)

`union(R1,R2)` or `R1 ∪ R2`\
Construct the union of two regions, returning a value of subtype `RegionUnion`. (Currently a stub for further development.)

There is also a parametric subtype `AbstractConnectedRegion{N}`, meant to represent a region of connectivity `N`. It has no associated method definitions.

## Generic types

### SimplyConnectedRegion

A `SimplyConnectedRegion` is a subtype of `AbstractConnectedRegion{1}`. The preferred construction is to call `interior(P)` or `exterior(P)` for an `AbstractClosedPath` or `AbstractClosedCurve` `P`. For bounded curves, these ignore the orientation of the boundary and select the bounded or unbounded regions, respectively. Otherwise, the points "to the left" of the boundary are considered the interior.

In addition to the methods of the [Abstract interface](@ref interface_regions), the type provides the following methods. (Here `R` is a `SimplyConnectedRegion`.)

`!(R)`\
Construct the region complementary to `R`. This is not a set complementation, as the boundary is not part of either region.

`isapprox(R1,R2)`\
Determine whether `R1` and `R2` represent the same region, regardless of boundary parameterization. This is done to within a tolerance that may be given to override the [global default](@ref global_defaults).

The `SimplyConnectedRegion` type is parameterized by the type of curve bounding it, to facilitate dispatch. Notably, there are definitions

```julia
AbstractDisk = SimplyConnectedRegion{T} where T<:Circle
AbstractHalfplane = SimplyConnectedRegion{T} where T<:Line
PolygonalRegion = SimplyConnectedRegion{Polygon}
```

There are also methods to facilitate construction of important common regions. For disks there are

- `disk(C::Circle)`
- `disk(center,radius)`
- `unitdisk`

For half-planes there are

- `halfplane(L::Line)`
- `upperhalfplane`
- `lowerhalfplane`
- `lefthalfplane`
- `righthalfplane`

### ConnectedRegion

The parameterized type `ConnectedRegion{N}` represents a region of connectivity `N`. You construct one by calling `ConnectedRegion{N}(outer,inner)`, where `outer` (if given) is an outer boundary, possibly unbounded, and `inner` is a vector of disconnected inner boundary components. Some rudimentary checking is done that a valid region of connectivity `N` has been specified, but it should not be considered rigorous.

The particular case of a doubly connected region can be constructed by `between(outer,inner)`, giving the two boundary components. The given orientation is ignored for any bounded component.

## [Specific subtypes](@id subtypes_regions)

### Annulus

An `Annulus` is the doubly connected region between two concentric circles. It is a subtype of `AbstractConnectedRegion{2}`. Construction is by `Annulus(outer,inner)`, where [Circle](@ref) values  are given explicitly, or by `Annulus(outrad,inrad,center=0)`, giving the radii and optionally the center.

## [Examples](@id examples_regions)
