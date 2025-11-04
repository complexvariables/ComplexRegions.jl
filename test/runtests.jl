using ComplexRegions, Statistics
CR = ComplexRegions
check(u::Number, v::Number, T::Type=real_type(float(u))) = isapprox(u, v, rtol=50CR.tolerance(T), atol=50CR.tolerance(T))
check(u::AbstractArray, v::AbstractArray) = all(check(u, v) for (u, v) in zip(u, v))
check(u::AbstractArray, v::Number) = all(check(u, v) for u in u)
using Test
@testset "Utilities" begin
    @test CR.scaleto(1im, 3im, [0.5, 0.75]) ≈ [2.0im, 2.5im]
    @test CR.scalefrom(1im, 3im, [2im, 1im, 1.5im]) ≈ [0.5, 0, 0.25]
    x = CR.realroots(3, 16, 1)
    @test all(@. abs(3x^2 + 16x + 1) < 1e-12)
    z = 2im .+ 3 * exp.(2im * pi * [0.8, 0.1, 0.25])
    @test CR.isccw(z...)
    @test CR.intadapt(exp, 0, 4, 1e-13) ≈ (exp(4) - 1)
    z = t -> cis(t)
    @test CR.fdtangent(z, 0) ≈ 1im
    @test CR.fdtangent(z, 1) ≈ 1im * exp(1im)
    @test CR.fdtangent(z, 0.2) ≈ 1im * cis(0.2)
end

@testset "Curves using $T" for T in (Float64, BigFloat)
    capprox(c1, c2) = all(check(point(c1, t), point(c2, t)) for t in (0:T(1000))/1000)
    f = t -> 2cospi(t) + 1im*sinpi(t)
    arc = Curve(f, -1, 0)
    @test arc isa Curve
    @test check(point(arc, 1//2), f(-1//2))
    @test check(point(conj(arc), 1//2), conj(f(-1//2)))
    @test check(point(inv(arc), 1//2), 1 / f(-1//2))
    arc = Curve(f)
    @test length(arc) == 1
    @test isfinite(arc)
    @test check(point(arc, 1//2), f(1//2))
    @test check(point(reverse(arc), 1//4), f(3//4))
    @test capprox(arc + 1im*T(1//3), Curve(t -> f(t) + 1im*T(1//3)))
    @test capprox(arc - 1im*T(1//3), Curve(t -> f(t) - 1im*T(1//3)))
    @test capprox(arc * 1im*T(1//3), Curve(t -> f(t) * 1im*T(1//3)))
    @test capprox(arc / 1im*T(1//3), Curve(t -> f(t) / 1im*T(1//3)))
    ellipse = ClosedCurve{T}(f, 0, 2)
    @test isfinite(ellipse)
    @test isfinite(reverse(ellipse))
    @test check(point(5 - 3im * ellipse, 1//8), 5 - 3im * f(T(1) / 4))
    @test check(point(conj(ellipse), 1//2), conj(f(1)))
    @test check(point(inv(ellipse), 1//2), 1 / f(1))
    @test check(angle(normal(ellipse, 3//4)), T(π) / 2)
    @test isapprox(arclength(ellipse), 9.6884482205476761984285031963918294, rtol=1e-10)
    @test winding(ellipse, 1 + 0.1im) == 1
    @test winding(ellipse)(1 - 2im) == 0
end

@testset "Circles in $T" for T in (Float64, BigFloat)
    c = Circle(Spherical{T}(1 - 1im), sqrt(T(2)))
    @test arclength(c) ≈ 2 * sqrt(T(2)) * T(pi)
    @test dist(-1 + 1im, c) ≈ sqrt(T(2))
    @test closest(1 + 4im, c) ≈ 1 + 1im * (sqrt(T(2)) - 1)
    @test isinside(1.5 - 1im, c) && isoutside(3//2 + 1im, reverse(c))
    @test !isinside(c)(2) && !isoutside(c)(-1//2*1im + 1)
    @test isinf(reflect(c.center, c))
    @test reflect(reflect(-1 + 2im, c), c) ≈ -1 + 2im
    @test all(mod(abs(arg(c, c(t)) - t), 1) < 1e-11 for t in 0:1//10:1)
    @test angle(unittangent(c, 1//8)) ≈ 3T(π)/4
    @test abs(tangent(c, T(1//8))) ≈ 2T(π) * sqrt(T(2))
    @test ispositive(c)
    @test length(c) == 1
    @test !ispositive(conj(c))
    @test !isapprox(c, Line(T(0), 1))
    τ = tangent(c, T(1//5))
    @test isapprox(τ, CR.fdtangent(c, T(1//5)), rtol=sqrt(eps(T)))

    c = Circle(1 - 1im, sqrt(2))
    @test point(c, 0.25) ≈ complex(1, sqrt(2) - 1)
    @test 2 / c isa Line
    c = Circle(1.0f0, -1im, 0)
    @test c.radius ≈ 1 / sqrt(2.0f0)
    c = Circle(1, -1im, false)
    @test all(mod(abs(arg(c, c(t)) - t), 1) < 1e-11 for t in 0:0.1:1)
    zz = point(5im - c / 3im, 0.23)
    @test abs(zz - (5im - c.center / 3im)) ≈ c.radius / 3
    @test Circle(1 + 3im, Polar(1 - 1im), 1.0) isa Line
end

@testset "Arcs in $T" for T in (Float64, BigFloat)
    a = Arc(cispi.([T(1) / 2, T(1) / 5, T(0)])...)
    zz = 1 / sqrt(T(2)) * (T(1) + 1im)
    @test check(point(a, T(1)/2) , zz)
    @test all(check(arg(a, a(t)), t) for t in (0:T(10))/10)
    @test check(dist(3im + T(1)/2 * 2im * cispi(T(1)/5), 3im + 2im * a), 2 * T(1)/2)
    @test isapprox(tangent(a, T(1)/5), CR.fdtangent(a, T(1)/5), rtol=sqrt(eps(T)))
    @test check(angle(unittangent(a, T(1)/2)), -T(π)/4)
    a = (a - 2) / 1im
    @test all(check(arg(a, a(t)), t) for t in (0:T(10))/10)
    @test inv(inv(a)) ≈ a

    b = Arc(-1im, 1im, -T(1))
    @test check(point(b, T(2) / 3), 1im)
    @test check(point(T(1)/10 - 3im * b, T(2) / 3) , T(1)/10 - 3im * 1im)
    @test check(closest(5im, b), 1im)
    @test check(closest(2 - 5im, b + 2) , 2 - 1im)
    @test check(dist(-2 - 4im, b), abs(-T(2) - 4im + 1im))
    @test check(closest(-2 - 4im, b), -T(1)*1im)
    @test all(check(arg(b, b(t)), t) for t in (0:T(10))/10)
    @test !isapprox(a, b)

    b = reverse(b)
    @test check(point(b, T(1) / 3), 1im)
    @test all(check(arg(b, b(t)), t) for t in (0:T(10))/10)
    @test check(angle(tangent(b, T(2) / 3)) , -T(π) / 2)

    b = Arc(Circle{T}(2im, 2), 1//3, 3//5)
    z = point(b, [0, 1])
    @test inv(b) ≈ Segment(1 ./ z...)
    @test check(arclength(b), T(π) * 12 / 5)
    @test conj(b) ≈ Arc(Circle{T}(-2im, 2), -1//3, -3//5)
end

@testset "Lines in $T" for T in (Float64, BigFloat)
    @test Line(1, 5) isa Line
    l = Line(T(1)*1im, direction=1 + 2im)
    @test slope(l) ≈ T(2)
    @test angle(l) ≈ atan(T(2))
    @test angle(conj(l)) ≈ -atan(T(2))
    @test angle(reverse(l)) ≈ atan(T(2)) - π
    @test isleft(2im, l) && !isleft(0, l)
    @test !isright(2im, l) && isright(0, l)
    @test isinf(arclength(l))
    @test !isfinite(l)
    @test ispositive(l)
    @test all(check(arg(l, point(l, t)), t) for t in (1:T(9))/10)
    @test isnothing(arg(l, 1 + 2im))
    @test check(abs(unittangent(l, T(1)/3)), 1)
    dz = point(l, T(3) / 5) - point(l, T(1) / 10)
    @test angle(dz) ≈ angle(1 + 2im)
    @test 1 / l isa Circle
    zz = point(5im + 1im*l / T(3), T(23)/100)
    z0 = T(5)im - l.base / 3im
    @test angle(zz - z0) ≈ angle(tangent(l, T(1)/2) / 3im)
    @test tangent(l, T(2)/5) ≈ CR.fdtangent(l, T(2)/5) rtol = sqrt(eps(T))
    z = l(T(3)/10) + 1im * sign(l.direction)
    @test dist(z, l) ≈ 1
    @test dist(z + 1 - 2im, (l + 1) - 2im) ≈ 1
    @test closest(z, l) ≈ l(T(3)/10)
    @test reflect(z, l) ≈ l(T(3)/10) - 1im * sign(l.direction)
    @test all(check(arg(l, l(t)), t) for t in (0:T(9))/10)
    @test convert(Line{Float32}, l) ≈ Line(1.0f0im, direction=1+2im)
end

@testset "Segments in $T" for T in (Float64, BigFloat)
    s = Segment(1, T(3) + 5im)
    @test isleft(-1, s) && !isleft(2, s)
    zz = 2 + T(5)*1im / 2
    @test point(2 - 3im * s, T(1)/2) ≈ 2 - 3im * zz
    @test closest(4 + 6im, s) ≈ 3 + 5im
    @test dist(-1, s) ≈ 2
    @test angle(tangent(s, 2//3)) ≈ angle(s(3//5) - s(1//10))
    @test tangent(s, T(3)/4) ≈ CR.fdtangent(s, T(3)/4) rtol=sqrt(eps(T))
    z = s(7//10) + 1im * sign(s(9//10) - s(7//10))
    @test closest(z, s) ≈ s(7//10)
    @test reflect(z, s) ≈ s(7//10) - (z - s(7//10))
    @test all(check(arg(s, s(t)), t) for t in (0:10)//10)
end

@testset "Rays in $T" for T in (Float64, BigFloat)
    s = Ray(Polar(2, 0), T(π) / 2)
    @test !isfinite(s)
    @test check(sign(s), T(1)*1im)
    @test isinf(arclength(s))
    @test isleft(-1im, s) && !isleft(-1im, reverse(s))
    @test check(real(s(23//100)), 2)
    @test imag(s(9//10)) > imag(s(7//10))
    @test check(convert(Complex, closest(T(5)*1im, s)), 2 + 5im)
    @test all(check(arg(s, s(t)), t) for t in (T(0):10)/10)
    @test check(2angle(tangent(s, 1//10)) , π)
    @test tangent(s, T(1)/10) ≈ convert(Complex, CR.fdtangent(s, T(1)/10)) rtol=sqrt(eps(T))
    @test check(2angle(tangent(reverse(s), 1)) , -T(π))
    z = T(3//7) + 2im
    @test s + z ≈ Ray(s.base + z, T(π) / 2)
    @test s - z ≈ Ray(s.base - z, T(π) / 2)
    @test s * z ≈ Ray(z * s.base, T(π) / 2 + angle(z))
    @test inv(s) ≈ Arc(T(1//2), T(1//4) - 1im/T(4), 0)

    s = Ray(Spherical(T(2)*1im), T(π), true)
    @test conj(s) ≈ Ray(-T(2)*1im, -T(π))
    @test !isright(3im, s) && isright(3im, reverse(s))
    @test isright(-3im, s) && !isright(-3im, reverse(s))
    @test check(sign(s), T(1))
    @test check(imag(s(1//2)) , 2)
    @test real(s(3//10)) < real(s(2//5))
    @test tangent(s, T(3)/5) ≈ CR.fdtangent(s, T(3)/5) rtol=sqrt(eps(T))
    @test !isleft(4, s) && isleft(-1 + 3im, s)
    @test check(convert(Complex, closest(-4 + 1im, s)) , -4 + 2im)
    @test check(convert(Complex, closest(6, s)) , 2im)
    @test all(check(arg(s, s(t)), t) for t in (T(0):10)/10)
end

@testset "Intersections in $T" for T in (Float64, BigFloat)
    u = T(1)/5 + 1im/T(2)
    z = intersect(Circle(0, 1), Circle(u, 3//2))
    @test all(@. abs(z - 0) ≈ 1)
    @test all(@. abs(z - u) ≈ T(3)/2)
    z = intersect(Circle(T(0), 1), Circle(T(1) + 1im, 3//2))
    @test all(@. abs(z - 0) ≈ 1)
    @test all(@. abs(z - (1 + 1im)) ≈ T(3)/2)
    @test isempty(intersect(Circle(0, 1), Circle(u, T(1)/10)))
    @test isempty(intersect(Circle(0, 1), Circle(u, 6)))

    z = intersect(Line(T(1), direction=1im), Line(-T(1), direction=1 + 1im))
    @test any(@. z ≈ 1 + 2im)
    @test isempty(intersect(Line(T(1), direction=1im), Line(-T(2), direction=1im)))
    l = Line(T(2), direction=3 + 1im)
    @test intersect(l, l + 100CR.tolerance(T)) ≈ l
    @test isempty(intersect(Segment(1im*T(1), 1im*T(3)), Segment(1im*T(-2), T(0))))
    @test intersect(Segment(1im*T(1), 1im*T(3)), Segment(1im*T(2), T(0))) ≈ Segment(1im*T(1), 1im*T(2))
    @test intersect(Segment{T}(1im, 3im), Line{T}(0, 1im)) ≈ Segment{T}(1im, 3im)
    @test isempty(intersect(Segment{T}(1im, 3im), Line(T(2), direction=1im)))

    z = intersect(Segment(0, 1), Segment(T(2)/5 - 1im, T(7//10) + 2im))
    @test any(@. z ≈ 0.5)
    z = intersect(Segment(T(2) + 3im, 3 + 3im), Segment(-T(1) + 3im, 4 + 3im))
    @test z ≈ Segment(T(2) + 3im, 3 + 3im)
    z = intersect(2 - Segment(T(2) + 3im, 3 + 3im), 2 - Segment(3im, T(24//10) + 3im))
    @test z ≈ 2 - Segment(T(2) + 3im, T(24//10) + 3im)
    z = intersect(1im * Segment(T(2) + 3im, 3 + 3im), 1im * Segment(T(27//10) + 3im, T(34//10) + 3im))
    @test z ≈ 1im * Segment(T(27//10) + 3im, 3 + 3im)
    z = intersect(Segment(2 + 3im, 3 + 3im), Segment(T(27//10) + 4im, T(34//10) + 6im))
    @test isempty(z)

    p1, p2 = Ray(0, T(pi) / 4), Segment(T(1//2) - 1im, T(1//2) + 2im)
    for z in (intersect(p1, p2), intersect(p2, p1))
        @test any(@. z ≈ (T(1) + 1im) / 2)
    end
    @test isempty(intersect(Ray(0, -3T(pi) / 4), Segment(T(1//2) - 1im, 0.5 + 2im)))
    @test isempty(intersect(Segment(T(1//2) - 1im, T(1//2) + 2im), Ray(0, -3T(pi) / 4)))
    z = intersect(Ray(0, T(pi) / 4), Ray(2 + 2im, T(pi) / 4))
    @test z ≈ Ray(2 + 2im, T(pi) / 4)
    z = intersect(Ray(0, T(pi) / 6), Ray(-1, -T(pi) / 2))
    @test isempty(z)

    p1, p2 = 3 + Circle(T(0), 1), 3 + Line(T(1//2), T(1//2) + 3im)
    for z in (intersect(p1, p2), intersect(p2, p1))
        @test any(@. z ≈ (T(7)/2 + sqrt(T(3)) / 2 * 1im) || z[1] ≈ (T(7//2) - sqrt(T(3)) / 2 * 1im))
    end

    p1, p2 = Circle(-2im, 2), Line(3im, 3)
    for z in (intersect(p1, p2), intersect(p2, p1))
        @test isempty(z)
    end

    p1, p2 = 3 + Ray(T(1//2) + T(1//10)*1im, T(pi) / 2, true), 3 + Circle(0, T(1))
    for z in (intersect(p1, p2), intersect(p2, p1))
        @test any(@. z ≈ (T(7//2) + sqrt(T(3)) / 2 * 1im)) && length(z) == 1
    end

    p1, p2 = 3 + Circle(T(0), 1), 3 + Segment(T(1//2), T(1//2) + 3im)
    for z in (intersect(p1, p2), intersect(p2, p1))
        @test all(@. z ≈ (T(7//2) + sqrt(T(3)) / 2 * 1im) || z ≈ (T(7//2) - sqrt(T(3)) / 2 * 1im))
    end

    p1, p2 = Circle(-2im, T(2)), Segment(-1im, T(1//2) - 1im)
    for z in (intersect(p1, p2), intersect(p2, p1))
        @test isempty(z)
    end
end

@testset "Paths in $T" for T in (Float64, BigFloat)
    S = Segment{T}(1, 1im)
    A = Arc(1im, -T(1) + 0.5im, -1)
    P = Path([S, A, -S])
    @test check(point(P, [0, 1, 1.5, 2.5, 3]), [S(0), S(1), A(0.5), -S(0.5), -S(1)])
    Q = 1 - 3im * P
    @test check(Q(1.5), 1 - 3im * A(0.5))

    z = (T(2) + 2im) / 5
    p = Path([Arc(-1, -z, -1im), Arc(-1im, conj(z), 1), Arc(1, z, 1im)])
    θ = angles(p)
    @test θ[2] ≈ 0.78121408739537 rtol=1e-13
    @test check(θ[2], θ[3])
    p = ClosedPath([curves(p); Arc(1im, -conj(z), -1)])
    θ = angles(p)
    @test check(θ[1:3], θ[2:4])

    P = Path([S, 1im * S, -S])
    @test check(vertex(P, 3), -1)
    @test check(vertex(P, 4), -point(S, 1))
    @test_throws BoundsError vertex(P, 5)
    @test check(arclength(P), 3 * sqrt(T(2)))
    @test length(vertices(P - 3im)) == 4
    rP = reverse(P)
    @test all(check(point(rP, t), point(P, 3 - t)) for t in range(T(0), T(3), 17))
    @test !isreal(P)
    @test isa(reverse(P), Path)
    τ = unittangent(P, T(3)/2)
    @test check(τ, 1im*unittangent(S, T(1)/2))
    @test (P * 3im) / 3im ≈ P
    Pi = inv(P)
    @test all(dist(1 / point(Pi, t), P) < CR.tolerance(T) for t in range(T(0), T(3), ))
    @test check(dist(T(1//3) + 1im, P), dist(T(1//3) + 1im, P[1]))
    @test check(closest(T(1//3) + 1im, P), closest(T(1//3) + 1im, P[1]))
    @test Path([convert(Segment{Float32}, S), 1im * S, -S]) ≈ P

    P = ClosedPath([S, 1im * S, -S, -1im * S])
    @test check(vertex(P, 3), -1)
    @test check(arclength(P), 4 * sqrt(T(2)))
    @test isa(reverse(P), ClosedPath)
    @test all(point(P, [0, 1, 5//4, 2.5, 3, 4]) .≈ [S(0), S(1), 1im * S(1//4), -S(0.5), -S(1), S(0)])
    @test winding(P, -1//3 - T(1//2) * 1im) == 1
    @test isinside(-1//3 - T(1//2) * 1im, P)
    @test winding(P, -1//3 + T(3//2) * 1im) == 0
    @test !isinside(-1//3 + T(3//2) * 1im, P)
    @test Path([convert(Segment{Float32}, S), 1im * S, -S, -1im * S]) ≈ P

    C = CR.enclosing_circle([P, -2P])
    @test all(isinside(point(P, t), C) for t in range(T(0), T(4), 100))
    @test all(isinside(point(-2P, t), C) for t in range(T(0), T(4), 100))

    P = ClosedPath([S, 1im * S, Arc(-T(1),-2 - 1im, 1)])
    @test check(2 * angles(P)[2], π)
end

@testset "Polygons in $T" for T in (Float64, BigFloat)
    s = Segment{T}(2, 2im)
    p = Polygon([s, 1im * s, -s, -1im * s])
    @test ispositive(p)
    @test check(arclength(p), 8 * sqrt(T(2)))
    @test check(-4*angle(normal(p, 1.1 + length(p))), π)
    z = (-T(4) + 5im) / 10
    @test winding(p, z) == 1
    @test winding(reverse(p), z) == -1
    @test winding(p, -T(4) - 0.5im) == 0
    @test check(2angles(p), π)
    @test ispositive(p)
    @test !ispositive(reverse(p))
    @test all(side(p, k) ≈ sides(p)[k] for k in 1:4)

    p = Polygon([T(4), 4 + 3im, 3im, -2im, 6 - 2im, 6])
    @test arclength(p) ≈ (3 + 4 + 5 + 6 + 2 + 2)
    @test tangent(p, 2.3 - length(p)) ≈ -5im
    @test winding(p, 5 - 1im) == 1
    @test winding(p, -1) == 0
    @test check(sum(angles(p)), 4T(π))

    p = CircularPolygon([Arc(T(1), 2 + 1im, 1im), Segment{T}(1im, -1), Arc(-T(1), -0.5im, -1im), Segment{T}(-1im, 1)])
    @test all(winding(p, z) == 1 for z in [1 + 0.5im, 1.7 + 1im, 0, -1 + 0.05 * cispi(1 / 5), -1im + 0.05 * cispi(0.3)])
    @test all(winding(p, z) == 0 for z in [-0.999im, 0.001 - 1im, -0.999, -1.001, 1.001, 1.999im, 1000im])
    @test check(arclength(p), sum(arclength(c) for c in curves(p)))
    @test all(curves(inv(p)) .≈ inv.(curves(p)))
    @test ispositive(p)

    r = Rectangle(T(2), [1, 3], T(π) / 2)
    @test arclength(r) ≈ 16
    @test all(abs.(imag(vertices(r))) .≈ 1)
    @test mean(real(vertices(r))) ≈ 2
    @test all(2angles(r) .≈ π )
    @test convert(Polygon, r) ≈ r.polygon
    @test point(r, 1//3) ≈ r.polygon(1//3)
    @test r(1//3) ≈ r.polygon(1//3)
    @test inv(r) ≈ inv(r.polygon)
    @test rectangle(vertices(r)) ≈ r

    for op in (+, -)
        z = T(2//3) + 2im
        @test op(r, z) ≈ op(Polygon(r), z)
    end

    for r in (rectangle(-2im, T(3) + 4im), rectangle([0, T(3)], [-2, 4]))
        @test arclength(r) ≈ 18
        @test all(extrema(r)[1] .≈ (0, 3))
        @test all(extrema(r)[2] .≈ (-2, 4))
    end

    @test all(length(n_gon(T, n)) == n for n in 3:7)
end

@testset "Unbounded polygons in $T" for T in (Float64, BigFloat)
    p = Polygon([T(5), 4 + 3im, 3im, -2im, 6 - 2im, (-T(π) / 2, 0)])
    a = angles(p) / π
    @test check(a[6], -0.5)
    @test check(sum(a .- 1), -2)

    p = Polygon([(T(π) / 2, T(π) / 2), T(5), 4 + 3im, 3im, -2im, 6 - 2im])
    a = angles(p) / π
    @test abs(a[1]) < CR.tolerance(T)
    @test check(sum(a .- 1), -2)
    @test all(winding(p, z) == 1 for z in [1 + 2im, 5 - 1im, 5.5 + 6im])
    @test all(winding(p, z) == 0 for z in [-3im, 3 + 5im, 5.5 - 6im])
    @test isfinite(truncate(p))

    p = Polygon([(-T(π) / 2, T(π) / 2), T(7), 4 + 3im, 3im, -2im, 6 - 2im])
    a = angles(p) / π
    @test check(a[1], -1)
    @test check(sum(a .- 1), -2)
    @test all(winding(p, z) == 1 for z in [4, 7 - 2im, 9])
    @test all(winding(p, z) == 0 for z in [4 + 4im, 6 + 2im, 4 - 3im])

    p = Polygon([4 + 3im, T(7), (0, 0), 6 - 2im, -2im, 3im])
    a = angles(p) / π
    @test check(a[3], -2)
    @test check(sum(a .- 1), -2)
end

@testset "Discretization in $T" for T in (Float64, BigFloat)
    p = Polygon([T(4), 4 + 3im, 3im, -2im, 6 - 2im, 6])
    t, z = discretize(p, 200)
    @test eltype(t) == T
    @test real_type(first(z)) == T
    @test length(t) == 200
    @test all(abs.(z .- p(t)) .< CR.tolerance(T))
    dz = abs.(diff(z))
    @test mean(dz) ≈   0.1100035321168 rtol=max(1e-12,CR.tolerance(T))
    @test median(dz) ≈ 0.1092971087529 rtol=max(1e-12,CR.tolerance(T))

    c = Circle{T}(1im, 2)
    t, z = discretize(c, 200)
    @test length(t) == 200
    @test eltype(t) == T
    @test real_type(first(z)) == T
    @test z ≈ c.(t) atol = CR.tolerance(T)
    @test all(abs.(z .- 1im) .≈ 2)

    s = Segment(T(-4), -2)
    t, z = discretize(s, 200)
    @test length(t) == 201
    @test eltype(t) == T
    @test eltype(z) == T
    @test z ≈ s.(t) atol = CR.tolerance(T)

    Z = discretize(interior(p), 100)
    @test size(Z) == (100, 100)
    @test count(isnan, Z) == 1981
    Z = discretize(interior(p), 100; limits=[-1, 1, -2, 2])
    @test all(isnan(z) || (abs(imag(z)) <= 2 && abs(real(z)) <= 1) for z in Z)
    Z = discretize(exterior(p), 120, limits=[-5, 10, -8, 7])
    @test size(Z) == (120, 120)
    @test count(isnan, Z) == 1536
    @test all(isnan(z) || (-5 <= real(z) <= 10 && -8 <= imag(z) <= 7) for z in Z)
    Z = discretize(exterior(p), 120)
    @test size(Z) == (120, 120)
    @test count(isnan, Z) == 1989
end

@testset "Möbius in $T" for T in (Float64, BigFloat)
    z = [T(0), 1, 2 + 2im]
    w = [Inf, -1im, -T(1)]
    for order = 1:3
        circshift!(z, 1)
        circshift!(w, 1)
        f = Möbius(z, w)
        @test check(f.(z), w)
        g = inv(f)
        @test check(g.(w), z)
    end
    c = Circle(z...)
    l = Line(w[2], w[3])
    @test f(c) ≈ l && g(l) ≈ c
    d = exterior(c)
    h = halfplane(l)
    @test f(d) ≈ !h && g(h) ≈ !d
    @test f(!d) ≈ h && g(!h) ≈ d
    f = Mobius(c, l)
    @test f(c) ≈ l
    a = Arc(c, 1//4, 1//2)
    @test f(a) ≈ Segment(f(c(1//4)), f(c(3//4)))
    f = Möbius([T(2), 3 + im, -5], z)
    g = Möbius(w, [T(2), 3 + im, -5])
    h = f ∘ g
    @test check(h.(w), z)
    f = Möbius(T(1), 2, 3, 4)
    g = Möbius([1 2; 3 T(4)])
    @test all(check(f(z), g(z)) for z in [T(0), 1, 2 + 2im])
end

@testset "Regions in $T" for T in (Float64, BigFloat)
    C1 = Circle{T}(1im, 1)
    C2 = Circle{T}(-1, 1//3)
    E = ExteriorRegion([C1, C2])
    for s in (-1//3, 2.2im, -3im), op in (+, -, *)
        @test op(s, -1im) ∈ op(s, E)
        @test op(s, -2 + 1im) ∈ op(E, s)
    end
    @test 2.2im / (2-1im) ∈ E / (2-1im)
    @test length(innerboundary(E)) == 2
    @test isnothing(outerboundary(E))
    @test !isfinite(E)
    C3 = Circle{T}(0, 10)
    F = connected_region(C3, [C1, C2])
    for s in (-1//3, 2.2im, -3im), op in (+, -, *)
        @test op(s, -1im) ∈ op(s, F)
        @test op(s, -2 + 1im) ∈ op(F, s)
        @test !(op(s, -12) ∈ op(F, s))
    end
    @test length(innerboundary(F)) == 2
    @test outerboundary(F) ≈ C3
    @test between(C3, C2) isa InteriorConnectedRegion{2}
    @test between(C3, reverse(C2)) isa InteriorConnectedRegion{2}
    @test between(reverse(C3), reverse(C2)) isa InteriorConnectedRegion{2}
end

@testset "SC Regions in $T" for T in (Float64, BigFloat)
    c = Circle{T}(0, 1)
    @test disk(c) ≈ interior(c)
    @test isfinite(interior(c))
    @test !isfinite(exterior(c))
    for c in (c, reverse(c)), D in (interior(c), !exterior(c))
        @test in(0, D)
        @test !in(Inf, D)
        @test in(Inf, !D)
        @test !in(0, !D)
    end

    for c in (Circle{T}(0, 1),
        CircularPolygon([Arc(T(1), 2 + 1im, 1im), Segment{T}(1im, -1), Arc(-T(1), -0.5im, -1im), Segment{T}(-1im, 1)]),
        Polygon([T(4), 4 + 3im, 3im, -2im, 6 - 2im, 6]))
        z = CR.get_one_inside(c)
        @test in(z, interior(c))
        @test !in(z, exterior(c))
        @test z isa Complex{T}
    end

    r = intersect(interior(Circle{T}(0, 1)), exterior(Circle{T}(1, 1//2)))
    @test -T(1//2) + 1im/T(5) ∈ r
    @test !(T(9//10) + 1im/T(7) ∈ r)
    r = union(interior(Circle{T}(0, 1)), exterior(Circle{T}(1, 1//2)))
    @test -T(1//2) + 1im/T(5) ∈ r
    @test T(9//10) + 1im/T(7) ∈ r
    r = intersect(interior(Circle{T}(0, 1)), interior(Circle{T}(3//2, 1)))
    @test !(-T(1//2) + 1im/T(5) ∈ r)
    @test !(T(12//10) - 1im/T(7) ∈ r)
    @test T(9//10) + 1im/T(7) ∈ r
    r = union(interior(Circle{T}(0, 1)), interior(Circle{T}(3//2, 1)))
    @test -T(1//2) + 1im/T(5) ∈ r
    @test T(12//10) - 1im/T(7) ∈ r
    @test T(9//10) + 1im/T(7) ∈ r

    @test boundary(interior(Circle{T}(0, 1))) ≈ Circle{T}(0, 1)
    @test outerboundary(interior(Circle{T}(0, 1))) ≈ Circle{T}(0, 1)
    @test isnothing(innerboundary(interior(Circle{T}(0, 1))))
    @test boundary(exterior(Circle{T}(0, 1))) ≈ Circle{T}(0, 1)
    @test innerboundary(exterior(Circle{T}(0, 1))) ≈ Circle{T}(0, 1)
    @test isnothing(outerboundary(exterior(Circle{T}(0, 1))))

    ell = Polygon([T(-10)-10im, 10-10im, 10-9im, -9-9im, -9+10im, -10+10im])
    L = interior(ell)
    @test boundary(!L) ≈ ell
    # stress the get_one_inside function
    @test in(CR.get_one_inside(ell), L)
    @test !in(CR.get_one_inside(ell), exterior(ell))
    @test in(CR.get_one_inside(ell + T(19)/2), 2L - T(-19))

    el = Line{T}(0, 1)
    @test halfplane(el) ≈ interior(el)
    @test halfplane(reverse(el)) ≈ exterior(el)
    @test halfplane(T(0), 1) ≈ halfplane(el)
    for H in (interior(el), !exterior(el))
        @test in(1im, H)
        @test !in(-1im, H)
        @test in(-1im, !H)
        @test !in(1im, !H)
    end

    @test CR.quad(rectangle(T(-1),3+3im)) isa CR.InteriorSimplyConnectedRegion{T}
end

@testset "Annulus in $T" for T in (Float64, BigFloat)
    A = Annulus(1, T(1//3))
    @test modulus(A) ≈ T(1//3)
    @test all(!in(z, A) for z in (0, 0.1im, -1.1im))
    @test all(in(z, A) for z in (0.4im, -0.99im))
    zc = 21 - T(3)*1im / 11
    A = Annulus(1, T(1//3), zc)
    @test modulus(A) ≈ T(1//3)
    @test all(!in(z + zc, A) for z in (0, 0.1im, -1.1im))
    @test all(in(z + zc, A) for z in (0.4im, -0.99im))
    @test isfinite(A)
    C1 = Circle{T}(zc, 1)
    C2 = Circle{T}(zc, T(1//3))
    A = Annulus(C1, C2)
    @test outerboundary(A) ≈ C1
    @test innerboundary(A)[1] ≈ C2
end

@testset "Shapes" begin
    @test all(θ ≈ π/3 for θ in angles(Shapes.triangle))
    @test all(θ ≈ π/2 for θ in angles(Shapes.square))
    @test length(Shapes.star) == 10
    @test length(Shapes.cross) == 12
    @test all(length(Shapes.hypo(k)) == k for k in 3:6)
    @test all(sum(reim(point(Shapes.squircle, t)).^4) ≈ 1 for t in 0:0.05:1)
    @test length(Shapes.spiral(3, 1)) == 4
    @test arclength(Shapes.spiral(3, 1)) ≈ 172.17394523
end

@testset "Utilities" begin
    @test isempty(CR.realroots(1, 2, 2))
    @test length(CR.realroots(1, -4, 4)) == 1
    @test length(CR.realroots(1, -4, -4)) == 2
    v = [-BigFloat(10) + 2im, 4 + 4im, 3 - 5im]
    for x in [1, 2]
        bre, bim = CR.enclosing_box(v, x)
        @test all(bre[1] .<= real(v) .<= bre[2])
        @test all(bim[1] .<= imag(v) .<= bim[2])
        zc, r = CR.enclosing_circle(v, x)
        @test all(abs.(v .- zc) .<= r)
    end
    for s in (Shapes.ellipse(2, 1/3), Shapes.squircle, Shapes.triangle, Shapes.cross)
        z = CR.adaptpoints(t -> s(t), t -> tangent(s, t), 0, 1)
        @test maximum(abs, diff(z)) < 0.03
    end
end
