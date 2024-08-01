# ComplexRegions

```@setup examples
using ComplexRegions
```

This package provides types and methods that are useful for working with curves and regions in the (extended) complex plane.

**Note:** Starting in version 0.2, plot capabilities have been moved to the [`ComplexPlots`](https://complexvariables.github.io/ComplexPlots.jl/stable/) package. Some of these are still used here to illustrate the examples.

Most functionality is provided through Julia types. Per Julia conventions, these are all capitalized. You use these capitalized names to create values of the type; e.g., [Segment](@ref) and [Circle](@ref).

Other methods may create values of these types, but since they are not distinct types themselves, they are not capitalized. For example, the [`rectangle`](@ref) method creates a [Polygon](@ref).

The methods in this package should work not only with the built-in `Complex` type, but also with the `Polar` and `Spherical` types from the [`ComplexValues`](https://complexvariables.github.io/ComplexValues.jl/stable/) package, which it re-exports.

## Abstract vs concrete types

All abstract types have names starting with `Abstract`. You probably won't encounter them unless you want to extend the provided functionality.

Abstract types cannot be instantiated. They serve as supertypes that collect common-denominator functionality, much like an interface or abstract class does in some object-oriented languages. Only the concrete descendants of these abstract types can be instantiated.

For example, any `AbstractCurve` is supposed to implement functions for finding points, tangents, and normals along the curve. There is a generic concrete `Curve` type that does the minimum required.  Specific subtypes such as a [Ray](@ref) or [Arc](@ref) provide (considerable) additional specialized functionalities appropriate to the subtypes.

## Curve, Path, and Region

A **curve** is meant to be a smooth, non-self-intersecting curve in the extended complex plane. The generic [Curve](@ref) type requires you to specify an explicit parameterization that is not checked for smoothness or even continuity. It will use automatic differentiation to find a tangent, if no tangent function is supplied. Particular subtypes of curve are `Circle`, `Arc`, `Line`, `Ray`, and `Segment`.

A **path** is a piecewise-continuous complex-valued path. In practice a [Path](@ref) can be specified as a vector of curves. The path is checked for continuity at creation time. The most important provided specific path types are [Polygon](@ref) and [CircularPolygon](@ref).

Both curves and paths have **closed** variants. These are additionally checked at creation to ensure that the initial and final points are the same.

One atypical aspect of curves and paths, even "closed" ones, is that they lie in the extended or compactified complex plane and thus may be unbounded. For instance, a line in the plane may be interpreted as a circle on the Riemann sphere, and is thus a closed curve passing through infinity.

Finally, a **region** is an open region in the extended plane bounded by a closed curve or path.

## Tolerance

Boundaries and endpoints are not well-posed ideas in floating-point, since an arbitrarily small perturbation to a value can move a point on or off of them. Thus many concepts in the package such as intersection or continuity are checked only up to a small tolerance. This value defaults to a modest multiple of machine precision and can be overridden on a per-call basis.