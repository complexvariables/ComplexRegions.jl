# Möbius

```@setup examples
using ComplexRegions
```

A **Möbius** transformation (also called bilinear or fractional-linear transformation) is the ratio of two linear polynomials:

$$f(z)=\frac{az+b}{cz+d}.$$

Among other notable properties, they map circles and lines to other circles and lines.

For convenience of typing on some keyboards, `Mobius` is an alias for `Möbius` in the package.

The package defines a `Möbius` type that can be constructed in a variety of ways:

- `Möbius(a,b,c,d)`\
Specify the coefficients as in the formula above.
- `Möbius(A)`\
Specify the coefficients as the matrix $A=[a\;b;\;\,c\;d]$.
- `Möbius(z,w)`\
Construct the unique transformation that maps the three points `z[1],z[2],z[3]` to `w[1],w[2],w[3]`, respectively. Either vector of points may include `Inf`.
- `Möbius(C1,C2)`\
Construct a transformation that maps the [Line](@ref) or [Circle](@ref) `C1` to the Line or Circle `C2`.

## Methods

Suppose `f` is a value of type `Möbius`. Then `f(z)` evaluates the transformation at the number `z`. In addition, `f(C)`, where `C` is a Circle or Line, returns the Circle or Line that is the image of `C` under `f`. Similarly, `f(R)`, where `R` is an `AbstractDisk` or `AbstractHalfplane`, returns the appropriate type of image region. For example,

```@repl examples
f = Möbius(Line(-1,1),Circle(0,1))
f(upperhalfplane)
isapprox(ans,unitdisk)
```

Two other methods are defined:

- `inv(f)`\
Construct the inverse transformation.
- `f∘g` (type "\circ" followed by tab key)\
Construct the composed map, $z \mapsto f(g(z))$.
