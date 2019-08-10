# Paths

A **path** is a sequence of curves that compose a continuous complex-valued path.

## [Abstract interface](@id interface_paths)

Every `AbstractPath` type is expected to provide the following method. (Here `P` is any value of type `AbstractPath`.)

`curves(P)`\
Return an array of the curves constituting `P`.

In addition, default implementations are given for the following methods.

`curve(P,k::Integer)`\
Return the `k`th curve of `P`.

`vertices(P)` or `vertex(P,k::Integer)`\
Return an array of the vertices of `P` (the endpoints of the constituent curves), or return only the `k`th vertex.

`length(P)`\
Number of curves in the path.

`isfinite(P)`\
True if the path is bounded.

`point(P,t::Real)` or `point(P,t::AbstractVector)`\
Compute a point or vector of points on the path `P`. For a path of length $n$, $t$ should be in $[0,n]$. The fractional part of $t$ corresponds to the parameter value on curve $1+\lfloor{k}\rfloor$ of the path, except at $t=n$, which is the final endpoint of the $n$th curve. 

`tangent(P,t::Real)` or `unittangent(P,t::Real)`\
Compute the complex tangent, or unit tangent, at a parameter value. The parameter is interpreted as in the `point` method. Note that the tangent at a vertex is not well-defined. 

`normal(P,t::Real)`\
Compute a complex, leftward-pointing unit normal at a parameter value. The parameter is interpreted as in the `point` method. Note that the normal at a vertex is not well-defined.

`conj(P)`\
Construct the complex conjugate of the path.

`reverse(P)`\
Reverse the orientation of the path.

`+,-,*,/`\
Translate, rotate and scale a path. If a path is in the denominator of a division, then it is inverted through the origin.

`isapprox(P1,P2)`
Determine whether two values represent the same path, within a tolerance.

`dist(z,P)`, `closest(z,P)`\
Respectively, find the distance from a point to the path, or find the point on the path nearest to a given number.

There is also an `AbstractClosedPath` subtype. It modifies a few of the implementations above:

`vertices(P)`\
Only returns the unique vertices; i.e., does not duplicate the initial/final vertex.

`curve(P,k)`, `vertex(P,k)`\
Returns the `k`th curve of the path in a circular/modulo sense; e.g., `k=0` is the same as `k=length(P)`.

`point`, `tangent`, `unittangent`, `normal`\
Evaluate using a circular interpretation of the parameter. E.g., if the length of the path is $n$, then the interval $[0,1]$ is equivalent to $[n,n+1]$.

### Iterator interface

The `AbstractPath` type implements the `eltype`, `length`, `getindex`, and `iterate` methods of the standard Julia iterator interface. Therefore, clauses such as `for c in P` will iterate over the curves in `P` for loops, comprehensions, and generators.

## Generic types

There are generic implementations of both of the abstract types described above.

### Path

`Path` implements the `AbstractPath` type. A path is created by calling `Path(c)`, where `c` is a vector of curves subtyped from `AbstractCurve`. The constructor tests the endpoints of the given curves for continuity up to a selectable tolerance.

In addition to the methods listed in the [Abstract interface](@ref interface_paths) above, a `Path` offers the `arclength` method. It requires all of the curve types it contains to have implementations of this same method, which is true for the built-in specific types, but not the `Curve` type.

### ClosedPath

`ClosedPath` implements the `AbstractClosedPath` type. Its chief difference from the `Path` type is that the constructor also checks whether the initial and final points coincide (up to tolerance). As a subtype of `AbstractClosedPath`, this type also inherits the use of circular/modulo addressing for curve selection and parameterization.

## Specific subtypes

There are two specialized implementations of the `AbstractClosedPath` type: `CircularPolygon` and the subtype `Polygon`, which implement the `AbstractCircularPolygon` and `AbstractPolygon` types, respectively.

### CircularPolygon

A `CircularPolygon` is a closed path whose curve components are all of type `Arc`, `Ray`, and `Segment`. In contrast to the usual notion of a polygon, the path may be unbounded. Construct a value by calling `CircularPolygon(c)` with a vector or `AbstractPath` of curves of appropriate types; continuity and closure of the path are checked as necessary.

In addition to the usual methods for a `ClosedPath`, the following are implemented.

`side`\
Alias for `curve`.

`winding(z,P)`\
Compute the winding number of `P` relative to `z`. Each counterclockwise rotation about `z` counts +1, and each clockwise rotation counts -1. The result is unreliable for points lying on `P`, for which the quantity is ill-posed.

`truncate(P)`\
Replace pairs of rays that meet at infinity with two segments to new vertices along the rays, and an arc between the new vertices. This can be useful for plotting an unbounded path, or for some other computations. This is *not* a true clipping algorithm; in fact, any bounded polygon will be unchanged.

### Polygon

A `Polygon` is a closed path whose curve components are all of type `Ray` and `Segment`. In contrast to the usual notion of a polygon, the path may be unbounded. Construct a value by calling `Polygon(c)` with a vector or `AbstractPath` of curves of appropriate types; continuity and closure of the path are checked as necessary.

An alternative construction is to provide a vector of vertices. In place of an infinite vertex, you can supply a tuple of the angles of the two rays that meet there. See the [Examples](@ref examples_paths) below. 

In addition to the methods for the [Abstract interface](@ref interface_paths) and [CircularPolygon](@ref), the `Polygon` type offers:

`angles(P)`\
Compute a vector of the interior angles of the polygon. Angles at a finite vertex are in the interval $(0,2\pi]$, while angles at an infinite vertex are in $[-2\pi,0]$, representing the angle at the pole of the Riemann sphere. 

## [Examples](@id examples_paths)

For example,
the following produces a polygon with two infinite vertices, resembling an infinite channel with a step:

```@example 1
using ComplexRegions # hide
p = Polygon([0,-1im,(0,0),1im,(pi,pi)])
```

```@example 1
using Plots # hide
plot(p)
savefig("channel.svg"); nothing # hide
```

![](channel.svg)
