using RecipesBase,Colors

RecipesBase.debug(false)

@recipe function f(C::AbstractCurve,n=600)
    aspect_ratio --> 1.0
    point(C,LinRange(0,1,n+1))
end

@recipe function f(C::Union{Line,Segment})
    aspect_ratio --> 1.0
    point(C,[0.0,1.0])
end

@recipe function f(P::AbstractPath;vertices=false)
    delete!(plotattributes,:vertices) 
    aspect_ratio --> 1.0  
    # hard coded in lieu of understanding plot themes
    palette --> [RGB{Float64}(0.0,0.6056031611752248,0.978680117569607)]
    label --> ""
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
