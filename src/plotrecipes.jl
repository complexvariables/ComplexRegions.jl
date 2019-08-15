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
            group := 2
            markercolor --> :black
            markershape --> :circle 
            seriestype := :scatter
            ComplexRegions.vertices(P) 
        end
    end
end

@recipe function f(p::AbstractCircularPolygon;vertices=false)
    if isfinite(p)
        p.path 
    else
        C = ComplexRegions.enclosing_circle(p,8)
        q = truncate(p,C)
        z,R = C.center,C.radius/4
        xlims --> (real(z)-R,real(z)+R)
        ylims --> (imag(z)-R,imag(z)+R)
        vertices := vertices
        q.path
    end
end

@recipe function f(R::SimplyConnectedRegion)
#    axis = plotattributes[:subplot][isvertical(plotattributes) ? :xaxis : :yaxis]
#    println(axis)
    seriestype := :shape
#    println(plotattributes[:plot_object])
    C = R.boundary  # could be curve or path
    if C isa Line 
        # need to fake with a polygon 
        θ = angle(C) 
        Polygon([C(0.5),(θ,θ+π)])
    elseif C isa AbstractCurve 
        ClosedPath(C) 
    else
        C 
    end
end
