# Curves

A **curve** is meant to be a smooth, non-self-intersecting curve in the extended complex plane. The curve parameter is a real number in the range $[0,1]$ that is used to traverse the curve.

All curve types are parameterized by a base floating-point type, such as Float64, that is the numeric type of the curve parameter as well as (possibly complexified) its points.

## [Abstract interface](@id interface_curves)

### AbstractCurve

Every realization of `AbstractCurve` is expected to implement the following methods. (Here `C` represents a value of type `AbstractCurve` and `z` is a number.)

| Method | Description |
|:-----|:-----|
| `point(C, t::Real)` | Complex point on `C` at parameter value `t` in [0,1].|
| `tangent(C, t::Real)` | Complex tangent to `C` at `t`.|
| `reverse(C)` | Reverse the direction of traversal.|
| `isfinite(C)` | True if the curve does not pass through infinity.|
| `conj(C)` | Complex conjugate of the curve. |
| `C+z` | Translate of the curve by `z`.|
| `-C` | Negate the curve.|
| `C*z` | Multiply the curve `C` by complex number `z`; i.e., scale and rotate it about the origin.|
| `inv(C)` | Invert the curve pointwise.|

There are also default implementations of the following methods:

| Method | Description |
|:-----|:-----|
| `point(C, t::AbstractArray{T<:Real})`| Vectorization of the `point` method. |
| `z+C, C-z, z-C, z*C, C/z, z/C` | Translate/rotate/scale by a complex value.|
| `unittangent(C, t::Real)`| Normalized tangent to `C` at `t`.|
| `normal(C, t::Real)`| Unit (leftward) normal to `C` at `t`.|
| `arclength(C)`| Arc length of `C`.|

### AbstractClosedCurve

The `AbstractClosedCurve` subtype is used to signify that the starting and ending points of the curve are (approximately) identical. In addition to the methods of `AbstractCurve`, it provides the following:

| Method | Description |
|:-----|:-----|
| `winding(C, z::Number)` | Winding number of `C` about `z`. |
| `isinside(z::Number, C)` | Detect whether `z` lies inside the curve. |
| `isoutside(z::Number, C)` | Detect whether `z` lies outside the curve. |

## [Generic types](@id generic_curves)

### Curve

A `Curve` represents an implementation of `AbstractCurve` that requires only an explicit parameterization of the curve. Given the (bounded) complex-valued function $f$ defined on $[0,1]$, then `Curve(f)` represents the curve $z=f(t)$. If $f$ is defined on $[a,b]$ instead, then `Curve(f,a,b)` is appropriate, but all future work with the curve object uses the standard interval $[0,1]$ for the parameter. All `Curve` values are expected to be finite; i.e., `isfinite(C)` will always be true.

By default, a tangent to the curve is computed when needed using automatic differentiation. If a function `df` is available for the complex-valued tangent $z'(t)$, it can be supplied via `Curve(f, df)`.

### ClosedCurve

A `ClosedCurve` implements `AbstractClosedCurve` and is similar to a `Curve`, but at construction the parameterization is checked for $f(0) \approx f(1)$, up to a tolerance.

## [Specific subtypes](@id subtypes_curves)

The following particular types of curves are provided. The default floating-point type is `Float64`, but you can specify another type explicitly, as in `Line{Float32}(0, 1im)` or `Segment(BigFloat(0), 1)`.

In addition to the minimal methods set by the `AbstractCurve` definition above, each of these types provides the following methods. (`C` is a value of one of these types, and `z` is a number.)

| Method | Description |
|:-----|:-----|
| `arg(C, z)`| Curve parameter value of a given point on the curve. |
| `isapprox(C1, C2)`| Determine whether two values represent the same curve.  |
| `isleft(z, C)`, `isright(z, C)`| Determine whether a point lies "to the left" or "to the right" of a line, ray, or segment in its given orientation. |
| `dist(z, C)` | Distance from a point to the curve. |
| `closest(z, C)`| Point on the curve nearest to a given number. |

### Line

Like other curves, a line is parameterized over $[0,1]$, with `L(0)` and `L(1)` both being infinity. 

- `Line(a, b)` creates the line through the values $a$ and $b$.
- `Line(p, direction=s)` creates a line through the point `p` in the direction of the complex number `s`.
- `Line(p, angle=θ)` creates a line through the point `p` at the angle `θ`.

Use `reflect(z, L)` to find the reflection of a point `z` across line `L`.

### Ray

- `Ray(z, θ)` constructs a ray starting at `z` and extending to infinity at the angle `θ`.
- `Ray(z, θ, true)` constructs a ray starting at infinity and extending to `z` at the angle `θ`.

### Segment

- `Segment(a, b)` constructs the line segment from `a` to `b`.

### Circle

- `Circle(z, r)` constructs a circle centered at `z` with radius `r`, oriented counterclockwise.
- `Circle(z, r, false)` constructs the circle with clockwise orientation.
- `Circle(a, b, c)` constructs the circle through the points `a`, `b`, and `c`. The ordering of the points determines the orientation of the circle. If the points are collinear, a `Line` is returned instead.

Use `reflect(z, C)` to reflect a point `z` through the circle `C`.

### Arc

- `Arc(a, b, c)` constructs the circular arc through the given three points. If the points are collinear, a `Segment` is returned.
- `Arc(C, start, Δ)` constructs an arc from a `Circle` `C`, starting at the given `start` value and extending an amount `Δ`. The values are expressed as fractions of a full rotation starting from the real axis.

## Examples

```@setup examples
using ComplexRegions
```

```@repl examples
ℓ = Line(1/2, 1/2+1im)    # line through 0.5 and 0.5+1i
c = 1 / ℓ    # a circle
winding(c, 1.5), winding(c, -1)
tangent(c, 0.75)
reflect(-1, c)
2c - 2
```
