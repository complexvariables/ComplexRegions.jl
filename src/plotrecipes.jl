using RecipesBase,Colors

RecipesBase.debug(false)

@recipe function f(C::AbstractCurve,n=600)
    aspect_ratio --> 1.0
    plotdata(C)
end

@recipe function f(P::AbstractPath;vertices=false,fillin=false)
    delete!(plotattributes,:vertices) 
    delete!(plotattributes,:fillin) 
    aspect_ratio --> 1.0  
    label --> ""
 
    data = vcat( [plotdata(c) for c in P]... )
    if fillin && isa(P,AbstractClosedPath)
        fillcolor --> :match 
        y = imag(vertex(P,1)+vertex(P,2))/2
        fillrange --> y
    end
    
    @series begin
        data
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
