using Documenter, ComplexRegions, Plots

makedocs(sitename="ComplexRegions.jl",
    format = Documenter.HTML(),
    authors = "Toby Driscoll",
    pages = [
        "Introduction" => "introduction.md",
        "Curves" => "curves.md",
        "Paths" => "paths.md",
        "Intersections" => "intersections.md",
        "Regions" => "regions.md",
        "API Reference" => "api.md"
		],
	modules = [ComplexRegions],
    doctest = true
    )

# deploydocs(
#     repo = "github.com/complexvariables/ComplexRegions.jl.git"
# #    versions = ["v#.#"],
# #    make = nothing
#     )
