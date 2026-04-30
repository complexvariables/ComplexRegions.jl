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
  curves.jl             # Abstract curve interface + Curve/ClosedCurve generics
  paths.jl              # Path/ClosedPath composing multiple curves
  regions.jl            # Abstract region types, set operations (∩, ∪, !)
  simplyconnected.jl    # Disk, half-plane, interior/exterior constructors
  polygons.jl           # Polygon, CircularPolygon, Rectangle, n_gon
  mobius.jl             # Möbius transformations
  lines.jl, circles.jl, arcs.jl, segments.jl, rays.jl  # Concrete curve types
  intersections.jl      # Pairwise curve intersection algorithms
  discretize.jl         # Arc-length equidistributed sampling
  utilities.jl          # intadapt, tolerance, enclosing_circle, etc.
ext/                    # Weak-dependency extensions for Plots, Makie, PythonCall
test/runtests.jl        # Comprehensive tests (Float64 and BigFloat)
```

## Type Hierarchy

```
AbstractCurve{T}  <: Function
├── AbstractClosedCurve{T}
│   ├── Line{T}, Circle{T}, ClosedCurve{T}
└── Segment{T}, Arc{T}, Ray{T}, Curve{T}

AbstractPath{T}
├── Path{T}
└── AbstractClosedPath{T}
    ├── ClosedPath{T}
    └── AbstractCircularPolygon{T}
        ├── CircularPolygon{T}
        └── AbstractPolygon{T} → Polygon{T}

AbstractRegion{T}
├── RegionIntersection, RegionUnion
└── AbstractConnectedRegion{N,T}   (N = connectivity)
    ├── ExteriorRegion, InteriorConnectedRegion
    ├── AbstractDisk, AbstractHalfplane
    └── AbstractSimplyConnectedRegion → Interior/ExteriorSimplyConnectedRegion
```

All types carry a real-type parameter `T` (typically `Float64` or `BigFloat`).

## Key Design Conventions

- **Curves are callable**: `C(t)` = `point(C, t)` with `t ∈ [0,1]`.
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
3. Default methods for `unittangent`, `normal`, `arclength`, `dist`, `closest`, `(C::MyType)(t)` are inherited.
4. Add algebraic operations (`+z`, `*z`, `inv`) following the pattern in existing types.
5. Add intersection methods in `intersections.jl` for pairwise combinations.

## Dependencies

Core: `CircularArrays`, `ComplexValues`, `Dierckx` (splines), `ForwardDiff` (autodiff for tangents of generic `Curve{T}`), `LinearAlgebra`, `StaticArrays`, `Statistics`, `Reexport`.

Visualization (weak): `Makie`, `Plots`, `PythonCall` — loaded only when the user imports them.

## Numerical Notes

- Adaptive integration (`intadapt`) uses Simpson's rule with recursive refinement.
- Generic `Curve{T}` tangents are computed via `ForwardDiff`; concrete types provide analytic tangents.
- `BigFloat` is tested alongside `Float64`; precision-sensitive code must use `tolerance(T)`.
- `discretize(p; ds)` samples by arc length; `refine_discretization` adapts based on curvature.
