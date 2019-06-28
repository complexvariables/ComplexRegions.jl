using RecipesBase

RecipesBase.debug(false)

@recipe function f(C::AbstractCurve,n=500)
    aspect_ratio --> 1.0
    point(C,LinRange(0,1,n+1))
end

@recipe function f(C::Union{Circle,Arc})
    aspect_ratio --> 1.0
    point(C,LinRange(0,1,601))
end

@recipe function f(C::Union{Line,Segment})
    aspect_ratio --> 1.0
    point(C,[0.0,1.0])
end

@recipe function f(P::AbstractPath;vertices=false)
    delete!(plotattributes,:vertices)    
    for c in P 
        @series begin 
            c
        end 
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
