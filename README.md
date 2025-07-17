# ComplexRegions

## Julia package for paths and regions in the complex plane 

[![][docs-stable-img]][docs-stable-url]
[![codecov](https://codecov.io/github/complexvariables/ComplexRegions.jl/graph/badge.svg?token=4RNN2G3NWY)](https://codecov.io/github/complexvariables/ComplexRegions.jl) 

[![DOI](https://zenodo.org/badge/210664766.svg)](https://zenodo.org/badge/latestdoi/210664766) 
[![DOI](https://joss.theoj.org/papers/10.21105/joss.01811/status.svg)](https://doi.org/10.21105/joss.01811)

This package provides types and methods that are useful for working with curves and regions in the complex plane.

Most functionality is provided through Julia types. Per Julia conventions, these are all capitalized. You use these capitalized names to create values of the type; e.g., [Segment](@ref) and [Circle](@ref).

Other methods may create values of these types, but since they are not distinct types themselves, they are not capitalized. For example, the [`rectangle`](@ref) method creates a [Polygon](@ref).

The methods in this package should work not only with the built-in `Complex` type, but also with the `Polar` and `Spherical` types from the [ComplexValues](https://complexvariables.github.io/ComplexValues.jl/stable/) package, which it re-exports.

Please see the [documentation](https://complexvariables.github.io/ComplexRegions.jl/stable/) for more details.

## Installation

The package can be installed with Julia's package manager:

```julia
julia> import Pkg
julia> Pkg.add("ComplexRegions")
```
or just use ```] add ComplexRegions``` at the usual command prompt.

## Project Status

Please open an [issue][issues-url] if you encounter any problems.

[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-latest-url]: https://complexvariables.github.io/ComplexRegions.jl/latest

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://complexvariables.github.io/ComplexRegions.jl/stable

[travis-img]: https://travis-ci.com/complexvariables/ComplexRegions.jl.svg?branch=master
[travis-url]: https://travis-ci.com/complexvariables/ComplexRegions.jl

[issues-url]: https://github.com/complexvariables/ComplexRegions.jl/issues
