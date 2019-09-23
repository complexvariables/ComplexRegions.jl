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

@recipe function f(p::AbstractCircularPolygon)
    if isfinite(p)
        p.path 
    else
        C = ComplexRegions.enclosing_circle(p,8)
        q = truncate(p,C)
        z,R = C.center,C.radius/3
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

@recipe function f(R::Union{ConnectedRegion,ExteriorRegion})
    p0 = outerboundary(R) 
    p1 = innerboundary(R)
    z1 = [plotdata(p) for p in p1]
    if isnothing(p0)
        zc,rho = enclosing_circle(vcat(z1...),8)
        p0 = Circle(zc,rho)
        r = 0.2*rho 
        xlims --> [real(zc)-r,real(zc)+r]
        ylims --> [imag(zc)-r,imag(zc)+r]    
    end
    z0 = plotdata(p0)
    # This is not fast, but I don't see a shortcut...
    # find pairwise distances between components
    comp = [z1...,z0]
    n = length(comp)
    index = Array{Tuple}(undef,n,n)
    dist = fill(Inf,n,n)
    for i in 1:n 
        for j in i+1:n
            ka,kb = argclosest(comp[i],comp[j])
            index[i,j] = (ka,kb)
            index[j,i] = (kb,ka)
            dist[i,j] = dist[j,i] = abs(comp[i][ka]-comp[j][kb])
        end
    end

    # find the hops between components  
    unused = trues(length(z1))
    curr = n
    path = []
    while sum(unused) > 1
        u = findall(unused)
        k = argmin(dist[curr,u])
        next = u[k]
        push!(path,(curr,next))
        unused[next] = false
        curr = next 
    end
    push!(path,(curr,findfirst(unused)))

    # accumulate into the last component 
    p = path[1] 
    idx = index[p...]
    data_in = z0[idx[1]] 
    data_out = [ z0[idx[1]:-1:1]; z0[end:-1:idx[1]] ]
    for k = 1:length(path)-1
        a = idx[2]
        zc = comp[p[2]]
        p = path[k+1]
        idx = index[p...]
        b = idx[1]
        if a > b
            data_in = [data_in;zc[a:-1:b]]
            data_out = [data_out;zc[a:end];zc[1:b]]
        else
            data_in = [data_in;zc[a:-1:1];zc[end:-1:b]]
            data_out = [data_out;zc[a:b]]
        end
    end
    a = idx[2]
    zc = comp[p[2]]
    data = [ data_in; zc[a:-1:1]; zc[end:-1:a]; data_out[end:-1:1] ]
 
    aspectratio --> 1
    @series begin
        seriestype := :shape
        linealpha := 0
       data
    end

    seriestype := :path
    linecolor --> :black 
    label := ""
    @series begin 
        z0
    end
    for z in z1
        @series begin
            z 
        end 
    end
end
