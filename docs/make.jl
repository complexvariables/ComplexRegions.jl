push!(LOAD_PATH,"../src")
import Pkg; Pkg.activate("..")
using Documenter, ComplexRegions

makedocs(sitename="ComplexRegions.jl",
    format = Documenter.HTML(),
    authors = "Toby Driscoll",
    pages = [
        "Home" => "index.md",
        "Introduction" => "introduction.md",
        "Curves" => "curves.md",
        "Paths" => "paths.md"
		],
	modules = [ComplexRegions],
    doctest = true
    )

# deploydocs(
#     repo = "github.com/complexvariables/ComplexRegions.jl.git"
# #    versions = ["v#.#"],
# #    make = nothing
#     )
