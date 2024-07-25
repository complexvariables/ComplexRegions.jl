# Use inverse linear interpolation to roughly equidistribute points on a curve.
function equidist!(t::AbstractVector, z::AbstractVector, p::AbstractCurveOrPath)
    s = [0; cumsum(abs.(diff(z)))]
    f = Spline1D(s, t, k=1)
    t .= f.(range(0, s[end], length=length(t)))
    @. z = p(t)
    return s[end]
end

function refine_discretization(p::AbstractCurveOrPath, lims::AbstractVector, ds::Real)
    t = collect(range(lims[1], lims[2], 20))
    idx = [1]
    while length(idx) > 0
        z = p.(t)
        S = equidist!(t, z, p)
        idx = findall(abs(z[i+1] - z[i]) > ds*S for i in 1:length(z)-1)
        t = [t; (t[idx] + t[idx.+1])/2]
        if length(t) > 50S / ds
            error("Too many points")
        end
        sort!(t)
    end
    if isclosed(p)
        t, z = t[1:end-1], z[1:end-1]
    end
    return t, z
end

"""
    discretize(p; ds=0.002, with_arg=false)
Discretize a path or curve, with points roughly equidistributed by arc length and chosen to be separated by no more than `ds` times the total arc length. All vertices are also included.

If `with_arg` is true, returns a tuple of vectors `t` and `z` such that `z[j]` is the point on the curve at parameter value `t[j]`. Otherwise, returns only `z`.
"""

function discretize(p::AbstractCurve{T}; ds=0.002, with_arg=false) where T
    lims = [T(0), T(1)]
    isinf(p(0)) && (lims[1] = T(1)/10)
    isinf(p(1)) && (lims[2] = T(9)/10)
    t, z = refine_discretization(p, lims, ds)
    return with_arg ? (t, z) : z
end

function discretize(p::AbstractPath{T}; ds=0.002, with_arg=false) where T
    t = T[]
    z = typeof(complex(p(T(1)/37)))[]
    for n in 1:length(p)
        tt, zz = discretize(p[n]; ds, with_arg=true)
        tt .+= n - 1
        if n > 1 && (T[1] == t[n-1][end])
            tt = tt[2:end]
            zz = zz[2:end]
        end
        append!(t, tt)
        append!(z, zz)
    end
    return with_arg ? (t, z) : z
end

"""
    discretize(p, n)
Discretize a path or curve at `n` points, roughly equidistributed by arc length. All vertices are also included.

Returns a tuple of vectors `t` and `z` such that `z[j]` is the point on the curve at parameter value `t[j]`.
"""
function discretize(p::AbstractClosedCurve{T}, n::Integer) where T
    t = isfinite(p) ? range(T(0), T(1), n+1) : range(T(1)/20, T(19)/20, n+1)
    t = collect(t)
    z = p.(t)
    equidist!(t, z, p)
    equidist!(t, z, p)
    return t[1:end-1], z[1:end-1]
end

function discretize(p::AbstractCurve{T}, n::Integer) where T
    t = isfinite(p) ? range(T(0), T(1), n+1) : range(T(1)/20, T(19)/20, n+1)
    t = collect(t)
    z = p.(t)
    equidist!(t, z, p)
    equidist!(t, z, p)
    return t, z
end

function discretize(p::AbstractPath{T}, n::Integer) where T
    m = length(p)
    isclosed(p) && (n += 1)
    t = isfinite(p) ? range(T(0), T(m), n) : range(T(1)/20, T(19m)/20, n)
    t = collect(t)
    z = p.(t)
    equidist!(t, z, p)
    equidist!(t, z, p)
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
    if isclosed(p)
        t, z = t[1:end-1], z[1:end-1]
    end
    return t, z
end

"""
    discretize(P::SimplyConnectedRegion, n=600)
Create an `n`×`n` grid of points on `P`. Points lying outside of `P` have a value of `NaN`.

If `P` is an exterior region, the points lie in a box a bit larger than the bounding box of `P`.

If keyword argument `limits` is specified, it must be a vector or tuple `(xmin, xmax, ymin, ymax)` specifying the limits of the grid.
"""
function discretize(
    P::InteriorSimplyConnectedRegion{T,S}, n=600;
    limits=nothing,
    ) where {T,S}
    @assert (isfinite(P) || !isnothing(limits)) "Unbounded region must have limits specified"
    # Get boundary points for determining interiority.
    _, z = discretize(boundary(P), 2n)
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
    T, z = discretize(boundary(P), 2n)

    if isnothing(limits)
        # Enlarge a bit to get an enclosing box.
        zc = mean(z)
        zzc = z .- zc
        r = max(maximum(real(zzc)), maximum(imag(zzc)))
        # zz = zc .+ 2.5*complex(r, r)
        # xlims, ylims = extrema(real(zz)), extrema(imag(zz))
        xlims = (real(zc) - 2r, real(zc) + 2r)
        ylims = (imag(zc) - 2r, imag(zc) + 2r)
    else
        xlims, ylims = Tuple(limits[1:2]), Tuple(limits[3:4])
    end

    # This function selects only outside points:
    point(x,y) = wind(complex(x,y), z) == 0 ? complex(x,y) : NaN
    return discretize(xlims, ylims, n, point)
end

function discretize(
    P::AbstractConnectedRegion, n=600;
    limits=nothing,
    )
    # Get boundary points for determining interiority.
    outer = outerboundary(P)
    zo = isnothing(outer) ? nothing : discretize(outer, 2n)[2]
    if isnothing(limits)
        xlims, ylims = extrema(real(zo)), extrema(imag(zo))
    else
        xlims, ylims = Tuple(limits[1:2]), Tuple(limits[3:4])
    end
    zi = map(p -> discretize(p, 2n)[2], innerboundary(P))

    # This function selects only inside points:
    function point(x, y)
        z = complex(x, y)
        if !isnothing(zo) && (wind(z, zo) == 0)
            return NaN
        else
            for c in zi
                if wind(z, c) != 0
                    return NaN
                end
            end
            return z
        end
    end

    return discretize(xlims, ylims, n, point)
end

# FIXME: kludgy
function discretize(E::ExteriorRegion, n::Integer=600)
    xlims = (Inf, -Inf)
    ylims = (Inf, -Inf)
    for c in innerboundary(E)
        z = discretize(c, ds=0.01)
        xx, yy = extrema(real(z)), extrema(imag(z))
        xlims = min(xlims[1], xx[1]), max(xlims[2], xx[2])
        ylims = min(ylims[1], yy[1]), max(ylims[2], yy[2])
    end
    r = max( xlims[2] - xlims[1], ylims[2] - ylims[1] ) / 2
    r *= 4//3
    xlims = mean(xlims) .+ (-r, r)
    ylims = mean(ylims) .+ (-r, r)
    return discretize(ConnectedRegion(nothing, E.inner), n, limits=(xlims..., ylims...))
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
