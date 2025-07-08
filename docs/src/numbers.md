# Number representations


```@setup examples
using ComplexRegions
```

The package is designed to work with complex numbers in different formats, including the native `Complex` type and the `Polar` and `Spherical` types from the [ComplexValues](https://github.com/complexvariables/ComplexValues.jl) package. The package re-exports these types, so you can use them directly.

```@example examples
seg = Segment(0, Polar(3, π/4))
point(seg, (0:4)/4)
```

```@example examples
Complex.(ans)
```

In addition, the package supports different floating-point types underlying any of the complex number types. This means you can use `BigFloats` or [`DoubleFloats`](https://github.com/JuliaMath/DoubleFloats.jl) to get higher precision.

```@example examples
using DoubleFloats
seg = Segment{DoubleFloat}(-1, 1)
seg(2//3)
```

Working in higher precision can be tricky: if you ever use a standard `Float64` value, it will set a ceiling on the precision of the result. For example, the following code will not give you the expected result:

```@example examples
seg(2/3)
```

As above, you can use `Rational` types to avoid premature floating-point conversion, or perform the conversions prior to calls:

```@example examples
seg(DoubleFloat(2) / 3)
```

Naturally, though, many powers of two can be converted correctly into higher precision without loss:

```@example examples
seg(1//8) - seg(0.125)
```

Curves, paths, and regions can be constructed with an explicit `AbstractFloat` type in braces. If the type is not given, the constructor will use the highest precision of its arguments, or default to `Float64` if the inputs are integers or rationals.

```@example examples
setprecision(100)
cir = Circle(BigFloat(1) / 5, 1//5)
```

Again, note the effect of premature float casting:

```@example examples
Circle(BigFloat(1) / 5, 1/5)
```
