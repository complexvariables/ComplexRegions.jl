module Shapes
using ComplexRegions
export star, triangle, square, cross, ellipse, hypo, circle, squircle, spiral

triangle = 1im * n_gon(3) / sqrt(3)
square = Polygon([-1-1im, 1-1im, 1+1im, -1+1im])

p = transpose(vertices(n_gon(5)))
star = Polygon( vec([p; cis(0.2π)*p/2]) )
star = 1im * star / vertex(star, 3)

function ellipse(a, b)
    r = t -> a*cospi(2t) + b*1im*sinpi(2t)
    dr = t -> 2π * (-a*sinpi(2t) + b*1im*cospi(2t))
    return ClosedCurve(r ,dr)
end
circle = ComplexRegions.Circle(0., 1.)

v = [-3-1im, -1-1im, -1-3im]
cross = Polygon([v; 1im*v; -v; -1im*v] / 3)

# Hypocycloids
function hypo(k::Integer)
    p = Curve[]
    z(t) = complex((k-1)*cospi(t) + cospi((k-1)*t), (k-1)*sinpi(t) - sinpi((k-1)*t))
    dz(t) = (2π/k) * complex(-(k-1)*sinpi(t) - sinpi((k-1)*t), (k-1)*cospi(t) - cospi((k-1)*t))
    for i in 0:k-1
        τ(t) = 2 * (i + t) / k
        push!(p, Curve(z∘τ, dz∘τ))
    end
    return ClosedPath(p)
end

function _squircle()
    ε = 1e-12   # regularize infinite tangents
    z(t) = complex(sqrt(cospi(t/2)), sqrt(sinpi(t/2)))
    dz(t) = (π/4) * complex(-1 ./ (ε + sqrt(sinpi(t/2))), 1 ./ (ε + sqrt(cospi(t/2))))
    s = Curve(z, dz)
    return ClosedPath([s, -conj(reverse(s)), -s, conj(reverse(s))])
end
squircle = _squircle()

"""
    spiral(n, w=0.5)
Create the boundary of a spiral region with `n` complete turns and width `w`.
"""
function spiral(n, w=0.5)
    z1 = t -> (1 + 2n*t) * cispi(2n*t)
    z1ʹ = t -> 2n * (1 + π * 1im * (1 + 2n*t)) * cispi(2n*t)
    z2 = t -> (1 + w + 2n*t) * cispi(2n*t)
    z2ʹ = t -> 2n * (1 + π * 1im * (1 + w + 2n*t)) * cispi(2n*t)

    return ClosedPath([
        Curve(z2, z2ʹ),
        Segment(z2(1), z1(1)),
        reverse(Curve(z1, z1ʹ)),
        Segment(z1(0), z2(0))
    ])
end

end
