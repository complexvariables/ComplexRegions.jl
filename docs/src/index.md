# ComplexValues

This package provides types and methods that are useful for working with curves and regions in the (extended) complex plane. 

# Examples

```@repl 1
using ComplexRegions
l = Line(1,1+1im)
1/l
c = Circle(-1im,2)
reflect(0,c)
intersect(l,c)
```
