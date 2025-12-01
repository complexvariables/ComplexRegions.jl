import Pkg; Pkg.update()
using DocumenterVitepress, Documenter
using ComplexRegions

ENV["GKSwstype"] = "100"
DocMeta.setdocmeta!(ComplexRegions, :DocTestSetup, :(using ComplexRegions); recursive=true)

makedocs(;
    modules = [ComplexRegions],
    repo = Remotes.GitHub("complexvariables", "ComplexRegions.jl"),
    authors = "Toby Driscoll <driscoll@udel.edu>, and contributors",
    sitename = "ComplexRegions.jl",
    format = DocumenterVitepress.MarkdownVitepress(
        repo = "https://github.com/complexvariables/ComplexRegions.jl",
    ),
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
        "Plotting in Makie" => "makie.md",
        "API Reference" => "api.md",
        "Using from Python" => "python.md",
    ],
)

DocumenterVitepress.deploydocs(;
    repo = "github.com/complexvariables/ComplexRegions.jl",
    target = "build", # this is where Vitepress stores its output
    devbranch = "master",
    branch = "gh-pages",
    push_preview = true,
)
