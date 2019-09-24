# ComplexRegions

```@setup examples
using ComplexRegions,Plots
default(linewidth=2,legend=:none)
```

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

One atypical aspect of curves and paths, even "closed" ones, is that they lie in the *extended* or compactified complex plane and thus may be unbounded. For instance, a line in the plane may be interpreted as a circle on the Riemann sphere, and is thus a "closed" curve passing through infinity.

A **region** is an open region in the extended plane bounded by a closed curve or path.

Some examples:

```@repl examples
ℓ = Line(1/2,1/2+1im)  # line through 0.5 and 0.5+1i
c = 1 / ℓ          # a circle
intersect(ℓ,c)
plot(ℓ);  plot!(c);
savefig("line_circle.svg"); nothing # hide
```

![line and circle](line_circle.svg)

```@repl examples
plot(Spherical(ℓ));
plot!(Spherical(c));
savefig("line_circle_sphere.svg"); nothing # hide
```

![line and circle on Riemann sphere](line_circle_sphere.svg)

```@repl examples
reflect(-1,c)       # reflection of a point through the circle
plot(interior(ℓ));   # plot a half-plane
savefig("halfplane.svg"); nothing # hide
```

![half-plan](halfplane.svg)

```@repl examples
h = n_gon(7)
plot(h);
for k in 1:7
	z = exp(k*2im*π/20)
	plot!(z*h - 0.5k - 0.1im*k^2)
end
savefig("heptagons.svg"); nothing # hide
```

![heptagons](heptagons.svg)

```@repl examples
p = Polygon([0,-1im,(0,0),1im,(pi,pi)])      # channel with a step
plot(interior(p));
savefig("channel.png"); nothing # hide
```

![channel with step](channel.png)

## Tolerance

Boundaries and endpoints are not well-posed ideas in floating-point, since an arbitrarily small perturbation to a value can move a point on or off of them. Thus many concepts in the package such as intersection or continuity are checked only up to a small tolerance. This value can be set on a per-call basis, or by using [global defaults](@ref global_defaults).

## [Global defaults](@id global_defaults)

For work at the REPL, it's convenient to be able to set an influential parameter just once rather than in multiple calls. This mechanism is provided via [`ComplexRegions.default`](@ref). You can see all the default parameters and values as follows:

```@repl examples
ComplexRegions.default()
```

Changing them is done with the same function:

```@repl examples
ComplexRegions.default(tol=1e-8)
```

Be advised that this type of "stateful" computing brings some subtle undesirable consequences. For example, if the global default `tol` is changed in a future release of the package, existing code could give different results when testing for interior points. If maximum reproducibility is a concern, you should develop the habit of setting all defaults yourself at the beginning of your code.
