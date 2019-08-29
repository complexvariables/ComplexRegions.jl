using RecipesBase,Colors

RecipesBase.debug(false)

@recipe function f(::Type{T},C::T) where T <: AbstractCurve
    aspect_ratio --> 1.0
    plotdata(C)
end

@recipe function f(P::AbstractPath;vertices=false)
    delete!(plotattributes,:vertices) 
    aspect_ratio --> 1.0  
 
    @series begin
        vcat( [plotdata(c) for c in P]... )
    end

    if vertices 
        @series begin 
            label := ""
            markercolor --> :black
            markershape --> :circle 
            seriestype := :scatter
            ComplexRegions.vertices(P) 
        end
    end
end

@recipe function f(R::AbstractConnectedRegion{2})
    p0 = outerboundary(R) 
    p1 = innerboundary(R)
    z0 = plotdata(p0)
    z1 = plotdata(p1)
    i1 = argmin( abs.(z1.-z0[1]) )
    @series begin
        aspectratio --> 1
        seriestype := :shape
        linealpha := 0
       [ z0[1];z1[i1:-1:1];z1[end:-1:i1];z0 ]
    end
    @series begin
        linecolor --> :black 
        p0 
    end
    @series begin
        linecolor --> :black 
        p1 
    end    
end

@recipe function f(p::AbstractCircularPolygon)
    if isfinite(p)
        p.path 
    else
        C = ComplexRegions.enclosing_circle(p,8)
        q = truncate(p,C)
        z,R = C.center,C.radius/4
        xlims --> (real(z)-R,real(z)+R)
        ylims --> (imag(z)-R,imag(z)+R)
        q.path
    end
end

@recipe function f(::Type{T},R::T) where T<:SimplyConnectedRegion
   if R isa ExteriorSimplyConnectedRegion 
        seriestype := :shapecomplement 
    else
        seriestype := :shape
    end

    C = R.boundary  # could be curve or path
    if C isa Line 
        # need to fake with a polygon 
        θ = angle(C) 
        Polygon([C(0.5),(θ,θ+π)])
#    elseif C isa AbstractCurve 
 #       ClosedPath(C) 
    else
        C 
    end
end

@recipe function f(R::ExteriorSimplyConnectedRegion)
    P = innerboundary(R)
    C = enclosing_circle(ClosedPath(P),8)
    zc = C.center
    r = 0.2*C.radius 
    xlims --> [real(zc)-r,real(zc)+r]
    ylims --> [imag(zc)-r,imag(zc)+r]
    between(C,P)
end
