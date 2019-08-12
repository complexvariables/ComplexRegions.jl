# ComplexRegions

This package provides types and methods that are useful for working with curves and regions in the (extended) complex plane.

Most functionality is provided through Julia types (roughly equivalent to classes in an object-oriented language). Per Julia conventions, these are all capitalized. You use these capitalized names to create values of the type; e.g., [Segment](@ref) and [Circle](@ref).

Other functions (methods, in Julia terms) may create values of these types, but since they are not distinct types themselves, they are not capitalized. For example, the [`rectangle`](@ref) method creates a [Polygon](@ref).

The methods in this package should work not only with the built-in `Complex` type, but also with the `Polar` and `Spherical` types from the [`ComplexValues`](https://complexvariables.github.io/ComplexValues.jl/stable/) package.

## Abstract types

All `abstract` types have names starting with `Abstract`. You probably won't encounter them unless you want to extend the provided functionality.

An abstract type cannot itself be instantiated as a value. They serve as supertypes that collect common-denominator functionality, much like an *interface* in other languages. For example, any `AbstractCurve` is supposed to provide functions for finding points, tangents, and normals along the curve. Specific subtypes such as a [Ray](@ref) or [Arc](@ref) provide additional specialized functionalities appropriate to the subtypes.

## Curve, Path, and Region

A **curve** is meant to be a smooth, non-self-intersecting curve in the extended complex plane. There is a generic [Curve](@ref) type that requires you to specify an explicit parameterization that is not checked for smoothness or even continuity. Implementations are given for more specific types of curve.

A **path** is a piecewise-continuous complex-valued path. In practice a [Path](@ref) can be specified as an array of curves. The path is checked for continuity at creation time. The most important provided specific path types are [Polygon](@ref) and [CircularPolygon](@ref).

Both curves and paths have **closed** variants. These are additionally checked that the initial and final points are the same.

A **region** is the open region in the extended plane bounded by a closed curve or path.

## Tolerance

Boundaries and endpoints are not well-posed ideas in floating-point, since an arbitrarily small perturbation to a value can move a point on or off of them. Thus many concepts in the package such as intersection or continuity are checked only up to a small tolerance. This value can be set on a per-call basis, or by using [global defaults](@ref global_defaults).

## [Global defaults](@id global_defaults)

For work at the REPL, it's convenient to be able to set an influential parameter just once rather than in multiple calls. This mechanism is provided via [`ComplexRegions.default`](@ref). You can see all the default parameters and values as follows:

```@repl 1
using ComplexRegions; # hide
ComplexRegions.default()
```

Changing them is done with the same function:

```@repl 1
ComplexRegions.default(tol=1e-8)
```

Be advised that this type of "stateful" computing brings some subtle undesirable consequences. For example, if the global default `tol` is changed in a future release of the package, existing code could give different results when testing for interior points. If maximum reproducibility is a concern, you should develop the habit of setting all defaults yourself at the beginning of your code.
