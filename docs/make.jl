import Pkg; Pkg.update("ComplexRegions")
using Documenter, ComplexRegions, Plots

DocMeta.setdocmeta!(ComplexRegions, :DocTestSetup, :(using ComplexRegions); recursive=true)
makedocs(sitename="ComplexRegions.jl",
    format = Documenter.HTML(),
    authors = "Toby Driscoll",
    pages = [
        "Introduction" => "index.md",
        "Curves" => "curves.md",
        "Paths" => "paths.md",
        "Polygons" => "polygons.md",
        "Intersections" => "intersections.md",
        "Regions" => "regions.md",
        "Möbius" => "mobius.md",
        "API Reference" => "api.md"
		],
	modules = [ComplexRegions],
    doctest = true
    )

deploydocs(
    repo = "github.com/complexvariables/ComplexRegions.jl.git"
#    versions = ["v#.#"],
#    make = nothing
    )
