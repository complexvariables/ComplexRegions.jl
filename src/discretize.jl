"""
    discretize(p::AbstractJordan, n=1200)
Discretize a closed path or curve at `n` points, roughly equidistributed by arc length. All vertices are also included.

Returns a tuple of vectors `t` and `z` such that `z[j]` is the point on the curve at parameter value `t[j]`.
"""
function discretize(p::AbstractClosedCurve, n=1000)
    t = range(0, 1, n+1)
    z = p.(t)
    # Use inverse linear interpolation to roughly equidistribute points.
    for _ in 1:2
        s = [0; cumsum(abs.(diff(z)))]
        f = Spline1D(s, t, k=1)
        t = f.(range(0, s[end], n+1))
        z = p.(t)
    end
    return t[1:end-1], z[1:end-1]
end

function discretize(p::AbstractClosedPath, n=min(1600, 500*length(p)))
    m = length(p)
    t = range(0, length(p), n+1)
    z = p.(t)
    # Use inverse linear interpolation to roughly equidistribute points.
    for _ in 1:2
        s = [0; cumsum(abs.(diff(z)))]
        f = Spline1D(s, t, k=1)
        t = f.(range(0, s[end], n+1))
        z = p.(t)
    end
    # Ensure that vertices are included.
    idx = 1
    for j in 1:m-1
        idx = findnext(>(j), t, idx)
        # replace the closest point
        v = p(j)
        if abs(v - z[idx-1]) < abs(z[idx] - v)
            idx -= 1
        end
        t[idx], z[idx] = j, v
    end
    return t[1:end-1], z[1:end-1]
end

"""
    discretize(P::SimplyConnectedRegion, n=600)
Create an `n`×`n` grid of points on `P`. Points lying outside of `P` have a value of `NaN`.

If `P` is an exterior region, the points lie in a box a bit larger than the bounding box of `P`.

If keyword argument `limits` is specified, it must be a vector or tuple `(xmin, xmax, ymin, ymax)` specifying the limits of the grid.
"""
function discretize(
    P::InteriorSimplyConnectedRegion, n=600;
    limits=nothing,
    )
    @assert isfinite(P) "Region must be bounded"
    # Get boundary points for determining interiority.
    z = discretize(boundary(P), 2n)[2]
    if isnothing(limits)
        xlims, ylims = extrema(real(z)), extrema(imag(z))
    else
        xlims, ylims = Tuple(limits[1:2]), Tuple(limits[3:4])
    end
    # This function selects only inside points:
    point(x,y) = wind(complex(x,y), z) != 0 ? complex(x,y) : NaN
    return discretize(xlims, ylims, n, point)
end

function discretize(
    P::ExteriorSimplyConnectedRegion, n=600;
    limits=nothing,
    )
    # Get boundary points for determining interiority.
    z = discretize(boundary(P), 2n)[2]

    if isnothing(limits)
        # Enlarge a bit to get an enclosing box.
        zc = mean(z)
        zz = zc .+ 2 * (z .- zc)
        xlims, ylims = extrema(real(zz)), extrema(imag(zz))
    else
        xlims, ylims = Tuple(limits[1:2]), Tuple(limits[3:4])
    end

    # This function selects only outside points:
    point(x,y) = wind(complex(x,y), z) == 0 ? complex(x,y) : NaN
    return discretize(xlims, ylims, n, point)
end

# Utility function for the main calls.
function discretize(xlims::NTuple{2}, ylims::NTuple{2}, n::Int, selector::Function)
    # selector(x,y) should return a complex number if (x,y) is inside P, and NaN otherwise.
    x = range(xlims..., length=n)
    y = range(ylims..., length=n)
    X = [x for x in x, y in y]
    Y = [y for x in x, y in y]
    Z = Matrix{ComplexF64}(undef, n, n)
    Threads.@threads for i in eachindex(Z)
        Z[i] = selector(X[i], Y[i])
    end
    return Z
end

# Fully discrete form of the winding number; faster than more precise versions
function wind(z0::Number, z::AbstractVector)
    n = length(z)
    q = 0
    for k in 1:n
        km = mod1(k-1, n)
        s = real(z[k]) - real(z0)
        sm = real(z[km]) - real(z0)
        if sign(s) != sign(sm)
            q -= sign(s - sm) * sign(imag(z[k]) - imag(z0))
        end
    end
    return q / 2
end
