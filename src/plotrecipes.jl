using RecipesBase,Colors

RecipesBase.debug(false)

@recipe function f(C::AbstractCurve,n=600)
    aspect_ratio --> 1.0
    label --> ""
    plotdata(C)
end

@recipe function f(P::AbstractPath;vertices=false)
    delete!(plotattributes,:vertices) 
    aspect_ratio --> 1.0  
    label --> ""
 
    @series begin
        vcat( [plotdata(c) for c in P]... )
    end

    if vertices 
        @series begin 
            markercolor --> :black
            markershape --> :circle 
            seriestype := :scatter
            vertex(P) 
        end
    end
end

@recipe function f(R::Region)
    fill --> true
    C = R.boundary
    C isa AbstractClosedCurve ? ClosedPath(C) : C
end
