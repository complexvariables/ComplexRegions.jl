import Pkg; Pkg.update("ComplexRegions")
using Documenter, ComplexRegions

ENV["GKSwstype"] = "100"
DocMeta.setdocmeta!(ComplexRegions, :DocTestSetup, :(using ComplexRegions); recursive=true)
makedocs(sitename="ComplexRegions.jl",
    format = Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://complexvariables.github.io/ComplexRegions.jl",
        edit_link="main",
        assets=String[],
    ),
    authors = "Toby Driscoll <driscoll@udel.edu>",
    repo="https://github.com/complexvariables/ComplexRegions.jl/blob/{commit}{path}#{line}",
    pages = [
        "Introduction" => "index.md",
        "Curves" => "curves.md",
        "Paths" => "paths.md",
        "Polygons" => "polygons.md",
        "Intersections" => "intersections.md",
        "Regions" => "regions.md",
        "Number types" => "numbers.md",
        "Möbius" => "mobius.md",
        "Shapes" => "shapes.md",
        "API Reference" => "api.md"
		],
	modules = [ComplexRegions],
    doctest = true
    )

deploydocs(
    repo = "github.com/complexvariables/ComplexRegions.jl",
    devbranch = "master",
#    versions = ["v#.#"],
#    make = nothing
    )
