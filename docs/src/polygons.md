# Polygons

There are two specialized implementations of the `AbstractClosedPath` type: `CircularPolygon` and the subtype `Polygon`, which implement the `AbstractCircularPolygon` and `AbstractPolygon` types, respectively.

## CircularPolygon

A `CircularPolygon` is a closed path whose curve components are all of type `Arc`, `Ray`, and `Segment`. In contrast to the usual notion of a polygon, the path may be unbounded. Construct a value by calling `CircularPolygon(c)` with a vector or `AbstractPath` of curves of appropriate types; continuity and closure of the path are checked as necessary.

In addition to the usual methods for a `ClosedPath`, the following are implemented:

- `side`\
Alias for `curve`.
- `winding(z,P)`\
Compute the winding number of `P` relative to `z`. Each counterclockwise rotation about `z` counts +1, and each clockwise rotation counts -1. The result is unreliable for points lying on `P`, for which the quantity is ill-posed.
- `truncate(P)`\
Replace pairs of rays that meet at infinity with two segments to new vertices along the rays, and an arc between the new vertices. This can be useful for plotting an unbounded path, or for some other computations. This is *not* a true clipping algorithm; in fact, any bounded polygon will be unchanged.

## Polygon

A `Polygon` is a closed path whose curve components are all of type `Ray` and `Segment`. In contrast to the usual notion of a polygon, the path may be unbounded. Construct a value by calling `Polygon(c)` with a vector or `AbstractPath` of curves of appropriate types; continuity and closure of the path are checked as necessary.

An alternative construction is to provide a vector of vertices. In place of an infinite vertex, you can supply a tuple of the angles of the two rays that meet there. See the [Examples](@ref examples_polygons) below.

In addition to the methods for the [Abstract interface](@ref interface_paths) and [CircularPolygon](@ref), the `Polygon` type offers

- `angles(P)`\
Compute a vector of the interior angles of the polygon. Angles at a finite vertex are in the interval $(0,2\pi]$, while angles at an infinite vertex are in $[-2\pi,0]$, representing the angle at the pole of the Riemann sphere.

Two additional special polygon constructors are defined:

- `rectangle(xlim,ylim)` or `rectangle(z1,z2)`\
Construct an axes-aligned rectangle given vectors of the real and imaginary limits, or two opposing complex corners.
- `n_gon(n)`\
Construct a regular n-gon with vertices on the unit circle.

# [Examples](@id examples_polygons)

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
