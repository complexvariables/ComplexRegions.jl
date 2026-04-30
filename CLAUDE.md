# ComplexRegions.jl

Julia package for curves and regions in the complex plane. Published in JOSS (v0.3.7). Used for conformal mapping research.

## Running Tests

```bash
julia --project=. -e "using Pkg; Pkg.test()"
# or within Julia REPL:
# ] test
```

## Project Layout

```
src/
  ComplexRegions.jl     # Main module — exports, type aliases, includes
  parameterized.jl      # AbstractParameterizedMap root (curves and paths share this)
  curves.jl             # Abstract curve interface + Curve/ClosedCurve generics
  paths.jl              # Path/ClosedPath composing multiple curves
  regions.jl            # Abstract region types, set operations (∩, ∪, !)
  simplyconnected.jl    # Disk, half-plane, interior/exterior constructors
  polygons.jl           # Polygon, CircularPolygon, Rectangle, n_gon
  shapes.jl             # Shapes module (square, star, ellipse, …)
  mobius.jl             # Möbius transformations
  lines.jl, circles.jl, arcs.jl, segments.jl, rays.jl  # Concrete curve types
  intersections.jl      # Pairwise curve intersection algorithms
  discretize.jl         # Arc-length equidistributed sampling
  utilities.jl          # intadapt, tolerance, enclosing_circle, etc.
ext/                    # Weak-dependency extensions for Plots, Makie, PythonCall
test/runtests.jl        # Comprehensive tests (Float64 and BigFloat)
```

## Type Hierarchy

Curves and paths share a common root that subtypes `Function`:

```
AbstractParameterizedMap{T} <: Function
├── AbstractCurve{T}
│   ├── AbstractClosedCurve{T}
│   │   ├── Line{T}, Circle{T}, ClosedCurve{T}
│   └── Segment{T}, Arc{T}, Ray{T}, Curve{T}
└── AbstractPath{T}
    ├── Path{T}
    └── AbstractClosedPath{T}
        ├── ClosedPath{T}
        └── AbstractCircularPolygon{T}
            ├── CircularPolygon{T}
            └── AbstractPolygon{T} → Polygon{T}
```

The `Jordan{T} = Union{AbstractClosedCurve{T}, AbstractClosedPath{T}}` alias is the dispatch surface for "any closed curve-or-path" — use `f(j::Jordan)` when smoothness/piecewise structure doesn't matter.

Regions form a separate tree with two abstract subdivisions of connectedness:

```
AbstractRegion{T}
├── RegionIntersection, RegionUnion
└── AbstractConnectedRegion{T}
    ├── InteriorRegion{T}
    ├── ExteriorRegion{T}
    ├── AbstractSimplyConnectedRegion{T}
    │   ├── InteriorSimplyConnectedRegion{T,S<:Jordan{T}}
    │   └── ExteriorSimplyConnectedRegion{T,S<:Jordan{T}}
    ├── AbstractDoublyConnectedRegion{T}
    │   ├── InteriorDoublyConnectedRegion{T,S<:Jordan{T}}
    │   ├── ExteriorDoublyConnectedRegion{T,S<:Jordan{T}}
    │   └── Annulus{T}
```

`connectivity(R::AbstractConnectedRegion)` returns the integer connectivity (1 for simply-connected, `length(innerboundary(R)) + (outer ? 1 : 0)` otherwise). Connectivity is *not* a type parameter — dispatch on the abstract types or check `connectivity(R) == n` at runtime when you need a specific value. Conformal-mapping algorithms specialized for the doubly-connected case can dispatch on `Annulus` directly or use the connectivity predicate.

Convenience aliases for common simply-connected regions:

- `Disk{T} = InteriorSimplyConnectedRegion{T, Circle{T}}`
- `Quad{T} = InteriorSimplyConnectedRegion{T, Rectangle{T}}`
- `Halfplane{T} = SimplyConnectedRegion{T, Line{T}}`
- `PolygonalRegion{T} = SimplyConnectedRegion{T, Polygon{T}}`

All types carry a real-type parameter `T` (typically `Float64` or `BigFloat`).

## Key Design Conventions

- **Curves are callable**: `C(t)` = `point(C, t)` with `t ∈ [0,1]`. Inherited from `AbstractParameterizedMap`.
- **Paths use integer intervals**: parameter `t ∈ [k, k+1]` maps to curve `k`.
- **Union type**: `AnyComplex{S} = Union{Complex{S}, Polar{S}, Spherical{S}}` — all three representations are accepted.
- **Tolerance**: `tolerance(T)` = `100 * eps(T)`. Use this (not hardcoded values) for floating-point comparisons.
- **Orientation**: Circle has a `ccw` flag; regions auto-normalize orientation in their constructors. Do not assume a curve is CCW without checking.
- **Algebraic operations**: `+`, `-`, `*`, `/`, `conj`, `inv` are defined on all curve types and return curves of compatible types.
- **Inside/outside**: use winding numbers via `winding(C, z)`, `isinside`, `isoutside`, or `z ∈ R` for regions.
- **Möbius**: `Möbius` struct stores `[a,b,c,d]` as a `SVector{4}` for `(az+b)/(cz+d)`. Maps circles/lines to circles/lines.

## Extending the Package

To add a new concrete curve type:
1. Define `struct MyType{T} <: AbstractCurve{T}` (or `AbstractClosedCurve`).
2. Implement: `point(C::MyType, t)`, `tangent(C::MyType, t)`, `reverse(C::MyType)`, `isfinite(C::MyType)`, `conj(C::MyType)`.
3. Default methods for `unittangent`, `unitnormal`, `arclength`, `dist`, `closest`, `(C::MyType)(t)` are inherited from `AbstractParameterizedMap`.
4. Add algebraic operations (`+z`, `*z`, `inv`) following the pattern in existing types.
5. Add intersection methods in `intersections.jl` for pairwise combinations.

To add a new region type, subtype the appropriate node:

- Simply-connected: `<: AbstractSimplyConnectedRegion{T}` and provide `boundary`, `in(z, R)`, `isfinite`.
- Multiply-connected: `<: AbstractMultiplyConnectedRegion{T}` and provide `outerboundary` (or `nothing`), `innerboundary` (a `Vector`), `in(z, R)`, `isfinite`.

## Dependencies

Core: `CircularArrays`, `ComplexValues`, `Dierckx` (splines), `ForwardDiff` (autodiff for tangents of generic `Curve{T}`), `LinearAlgebra`, `StaticArrays`, `Statistics`, `Reexport`.

Visualization (weak): `Makie`, `Plots`, `PythonCall` — loaded only when the user imports them.

## Numerical Notes

- Adaptive integration (`intadapt`) uses Simpson's rule with recursive refinement.
- Generic `Curve{T}` tangents are computed via `ForwardDiff`; concrete types provide analytic tangents.
- `BigFloat` is tested alongside `Float64`; precision-sensitive code must use `tolerance(T)`.
- `discretize(p; ds)` samples by arc length; `refine_discretization` adapts based on curvature.

## Known Gotchas

- Because `AbstractParameterizedMap <: Function`, broadcasting treats curves and paths as scalars (`Base.broadcastable(::Function) = Ref(x)`). So `tangent.(p::Polygon)` calls `tangent(p)` — a `MethodError` — rather than iterating sides. Use `tangent.(sides(p))` or a comprehension when you want per-side broadcast.
