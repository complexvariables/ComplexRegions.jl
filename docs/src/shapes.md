# Predefined shapes

```@setup examples
using ComplexRegions
```

The `Shapes` submodule provides a number of predefined shapes that can be used in constructing regions. 

```@example examples
using ComplexPlots, Plots
shapes = [
    Shapes.circle  Shapes.ellipse(2, 1) Shapes.squircle; 
    Shapes.square  Shapes.triangle      Shapes.cross;
    Shapes.hypo(3) Shapes.star          Shapes.spiral(2, 0.7)
    ]

fig = plot(resolution=(400, 400), layout=(3,3), showaxis=false, legend=false)
for i in 1:3, j in 1:3
    plot!(shapes[i, j], l=2, subplot=3i + j - 3, aspect_ratio=1, guide="", grid=false)
end
savefig("shapes.svg"); nothing # hide
```

![table of shapes](shapes.svg)
