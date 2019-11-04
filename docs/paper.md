---
title: 'ComplexRegions.jl: A Julia package for regions in the complex plane'
tags:
  - Julia
  - complex variables
  - conformal mapping
authors:
  - name: Tobin A. Driscoll
    orcid: 0000-0002-1490-2545
    affiliation: "1" # (Multiple affiliations must be quoted)
affiliations:
 - name: University of Delaware
   index: 1
date: 2 October 2019
bibliography: paper.bib
---

# Summary

Complex variables offer a unique and powerful way to represent planar regions and maps between them. Since the work of Riemann and others in the 19th century, complex variables have been used to derive quantities and properties of importance to a wide range of applications in science and engineering, in addition to their inherent mathematical importance. The advent of computing accelerated and extended this trend to applications in fluids, electromagnetics, queueing theory, computer graphics, and many other fields; e.g., [@DE92;@Dri94;@Gai79;@Gre90;@MBH92;@TW86;@Ver82;@Gu2004].

Software packages for working with complex variables, particularly conformal maps, have been in the public domain for decades [@BG87;@Crowdy2016;@Dri96;@HP83;@Hu95;@Marshall2007;@Tre80]. There are additional papers based on research codes that apparently have not made it into the public domain, e.g. [@DE93;@KT86;@OR89]. These packages tend to take the traditional form of library routines that require some degree of facility with custom data structures used to represent geometry and functions. Moreover, computational methods are mostly specialized to particular types or classes of regions, leading to parallel and duplicated functionality within a package, and incompatible representations between them.

The `ComplexRegions` package for Julia provides a software framework for the front end of computations over regions in the extended complex plane. It defines abstract types for curves, paths (i.e., piecewise smooth contours), and regions, defining minimal interfaces for them and providing default behavior based only on the interfaces. These are then implemented as concrete subtypes; e.g., a `Circle` is a subtype of `AbstractClosedCurve`, and a `Polygon` is a subtype of `AbstractClosedPath`. These types define data structures and methods to operate on them.

Julia's multiple dispatch facility enables some convenient uses of this basic framework. For example, the `Base.intersect` method is extended to have definitions for many possible pairings of curve and path arguments. In the future, a method for constructing conformal maps, say, could have one method for arguments of types `AbstractDisk` and `PolygonalRegion` that would call upon a Schwarz-Christoffel mapping code. The end user need not know what sort of mapping algorithm is needed for a particular problem, and the master method could easily be extended to numerous other contexts without the need for a "switchyard" portal that could become difficult to maintain.

The `ComplexRegions` package builds on the `ComplexValues` package that defines `Polar` and `Spherical` types for working with polar and Riemann sphere representations of complex numbers. It also provides recipes for plotting the major abstract types with the popular `Plots.jl` package. For example:
```
using ComplexRegions,Plots
c = Circle(0,1)
t = n_gon(3)
s = n_gon(4)
plot( ExteriorRegion([c,3+s,6+t]),leg=:none )
```
![A multiply connected exterior region](triple.pdf)

```
julia> li = Line(1/2,1/2+1im)  # line through 0.5 and 0.5+1i
Line{Complex{Float64}} in the complex plane:
   through (0.5 + 0.0im) parallel to (0.0 + 1.0im)

julia> c = 1 / li          # a circle
Circle{Complex{Float64}} in the complex plane:
   centered at (1.0 + 0.0im) with radius 1.0, negatively oriented

julia> intersect(li,c)
   2-element Array{Complex{Float64},1}:
    0.5 + 0.8660254037844386im
    0.5 - 0.8660254037844386im

julia> plot(Spherical(li),leg=:none);  plot!(Spherical(c))
```
![Line and circle on the Riemann sphere](line_circle.pdf)

[`DomainSets.jl`](https://github.com/JuliaApproximation/DomainSets.jl) and [`IntervalSets,jl`](https://github.com/JuliaMath/IntervalSets.jl), by D. Huybrechs, S. Olver, *et al.*, are a pair of packages with similar goals and interface, but oriented toward the representation of function domains in real spaces. Although `ComplexRegions` was not consciously patterned after them, there is enough similarity in approach to consider merging the functionality in the future.

One future use case for this package is to implement numerical methods for conformal maps. Another is for the fast solution of special cases of the Laplace equation that are especially amenable to solution using complex-variable methods [@GT19]. 

# References
