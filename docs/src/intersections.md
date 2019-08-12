# Intersections

Methods are provided to find intersections between all the [Specific subtypes](@ref subtypes_curves). In the generic cases where the intersections consist of zero or more points, a vector of results is returned. In special circumstances of partially or wholly overlapping curves, a Curve subtype is returned.

There are also methods for finding the intersections between curves and paths, and between two paths. These return set unions over the curves of the path(s), returning a vector with complex element type or, if some overlaps occurred, mixed types.

## [Examples](@id examples_intersections)


