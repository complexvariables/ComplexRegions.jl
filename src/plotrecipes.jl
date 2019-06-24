using RecipesBase

@recipe function f(z::Array{Polar{T}}) where T
    projection --> :polar
    Complex.(z)
end

@recipe function f(z::Array{Spherical{T}}) where T
    markersize --> 1
    x = [ cos(z.lat)*cos(z.lon) for z in z ]
    y = [ cos(z.lat)*sin(z.lon) for z in z ]
    z = [ sin(z.lat) for z in z ]
    x,y,z
end
