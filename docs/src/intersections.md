# Intersections

Methods are provided to find intersections between all the [Specific subtypes](@ref subtypes_curves). In the generic cases where the intersections consist of zero or more points, a vector of results is returned. In special circumstances of partially or wholly overlapping curves, a Curve subtype is returned.

There are also methods for finding the intersections between curves and paths, and between two paths. These return set unions over the curves of the path(s), returning a vector with complex element type or, if some overlaps occurred, mixed types.

## [Examples](@id examples_intersections)

```@setup examples
using ComplexRegions,Plots 
default(linewidth=2,markersize=3,legend=:none)
```

Circles and Arcs intersect at zero, one, or two points, or as an Arc. 

```@repl examples
c = Circle(0,1)
a = Arc(1+1im,0,1-1im)
intersect(c,a)
b = Arc(1im,1,-1im)
intersect(c,b)
ans ≈ b
```

Segments and Lines intersect at zero or one point, or as a Line or Segment. 

```@repl examples
l = Line(1im,1+1im)
s = Segment(-2,2+2im)
intersect(l,s)
intersect(l,s+2im)
intersect(s,Segment(-4-1im,1im))
```

Here is the intersection of a polygon with a segment.

```@repl examples
box = [1-1im,3-1im,3+1im];
plus = Polygon([box;1im*box;-box;-1im*box])
s = Segment(-3-2im,2+2im)
z = intersect(plus,s)
plot(plus);
plot!(s);
scatter!(Complex.(z));
savefig("plus_intersect.svg"); nothing # hide
```

![intersection of polygon and segment](plus_intersect.svg)
