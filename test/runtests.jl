using ComplexRegions
CR = ComplexRegions

using Test
@testset "Utilities" begin
	 @test(CR.scaleto(1im,3im,[0.5,0.75]) ≈ [2.0im,2.5im])
	 @test(CR.scalefrom(1im,3im,[2im,1im,1.5im]) ≈ [0.5,0,0.25])
	 x = CR.realroots(3,16,1)
	 @test( all(@. abs(3x^2 + 16x + 1)<1e-12 ) )
	 z = 2im .+ 3*exp.(2im*pi*[0.8,0.1,0.25])
	 @test( CR.isccw(z...) )
	 @test( CR.intadapt(exp,0,4,1e-13) ≈ (exp(4)-1)  )
	 z = t -> exp(1im*t)
	 @test( CR.fdtangent(z,0) ≈ 1im )
	 @test( CR.fdtangent(z,1) ≈ 1im*exp(1im) )
	 @test( CR.fdtangent(z,0.2) ≈ 1im*exp(0.2im) )
end

@testset "Curves" begin
	f = t -> 2*cos(t) + 1im*sin(t)
	@test( Curve(f,-1,1) isa Curve)
	@test( point(Curve(f,-1,1),0.5) ≈ f(0) )
	f = t -> 2*cos(t) + 3im*sin(t)
	c = ClosedCurve(f,0,2π)
	@test( point(5-3im*c,0.125) ≈ 5-3im*f(pi/4) )
	@test( angle(normal(c,0.75)) ≈ 0.5π )
end

@testset "Circles" begin
	c = Circle(Spherical(1-1im),sqrt(2))
	@test( arclength(c) ≈ 2*sqrt(2)*pi )
	@test( dist(-1+1im,c) ≈ sqrt(2) )
	@test( closest(1+4im,c) ≈ 1+1im*(sqrt(2)-1) )
	@test( isinside(1.5-1im,c) && isoutside(1.5+1im,reverse(c)) )
	@test( isinf(reflect(c.center,c)) )
	@test( reflect(reflect(-1+2im,c),c) ≈ -1+2im )
	@test( all( abs(arg(c,c(t))-t) < 1e-11 for t in 0.1:0.1:1 ) )
	@test( angle(unittangent(c,0.125)) ≈ 0.75π )
	@test( abs(tangent(c,0.125)) ≈ 2π*sqrt(2) )

	c = Circle(1-1im,sqrt(2))
	@test( point(c,.25) ≈ complex(1,sqrt(2)-1) )
	@test( 2/c isa Line )
	c = Circle(1f0,-1im,0)
	@test( c.radius ≈ 1/sqrt(2f0) )
	c = Circle(1,-1im,false)
	@test( all( abs(arg(c,c(t))-t) < 1e-11 for t in 0.1:0.1:1 ) )
	zz = point(5im - c/3im,.23)
	@test( abs(zz-(5im-c.center/3im)) ≈ c.radius/3 )
	@test( tangent(c,0.7) ≈ CR.fdtangent(c,0.7) )
	@test( Circle(1+3im,Polar(1-1im),1.0) isa Line )
end

@testset "Arcs" begin
	a = Arc(exp.(1im*[pi/2,pi/5,0])...)
	zz = 1/sqrt(2)*(1+1im)
	@test( point(a,0.5) ≈ zz )
	@test( all( abs(arg(a,a(t))-t) < 1e-11 for t in 0:0.1:1 ) )
	@test( dist(3im+.5*2im*exp(1im*pi/5),3im+2im*a) ≈ 2*0.5 )
	@test( tangent(a,0.2) ≈ CR.fdtangent(a,0.2) )
	@test( angle(unittangent(a,0.5)) ≈ -0.25π )
	a = (a-2)/1im 
	@test( all( abs(arg(a,a(t))-t) < 1e-11 for t in 0:0.1:1 ) )

	b = Arc(-1im,1im,-1)
	@test( point(b,2/3) ≈ 1im )
	@test( point(0.1 - 3im*b,2/3) ≈ 0.1 - 3im*1im )
	@test( closest(5im,b) ≈ 1im )
	@test( closest(2-5im,b+2) ≈ 2-1im )
	@test( all( abs(arg(b,b(t))-t) < 1e-11 for t in 0:0.1:1 ) )

	b = reverse(b) 
	@test( point(b,1/3) ≈ 1im )
	@test( all( abs(arg(b,b(t))-t) < 1e-11 for t in 0:0.1:1 ) )
	@test( angle(tangent(b,2/3)) ≈ -0.5π )
end 

@testset "Lines" begin
	@test( Line(1,5) isa Line )
	l = Line(1im,direction=1+2im)
	@test( isleft(2im,l) && !isleft(0,l) )
	dz = point(l,.6)-point(l,.1)
	@test( angle(dz) ≈ angle(1+2im) )
	@test( 1/l isa Circle )
	zz = point(5im - l/3im,.23)
	z0 = 5im - l.base/3im
	@test( angle(zz-z0) ≈ angle(tangent(l,0.5)/3im) )
	@test( tangent(l,0.2) ≈ CR.fdtangent(l,0.2) )
	z = l(0.3) + 1im*sign(l.direction)
	@test( dist(z,l) ≈ 1 )
	@test( closest(z,l) ≈ l(0.3) )
	@test( reflect(z,l) ≈ l(0.3) - 1im*sign(l.direction))
	@test( all( abs(arg(l,l(t))-t) < 1e-11 for t in 0:0.1:0.9 ) )
end

@testset "Segments" begin
	s = Segment(1,3.0+5.0im)
	@test( isleft(-1,s) && !isleft(2,s) )
	zz = 2 + 2.5im
	@test( point(2 - 3im*s,0.5) ≈ 2 - 3im*zz)
	@test( closest(4+6im,s) ≈ 3+5im )
	@test( dist(-1,s) ≈ 2 )
	@test( angle(tangent(s,2/3)) ≈ angle(s(0.6)-s(0.1)) )
	@test( tangent(s,0.75) ≈ CR.fdtangent(s,0.75) )
	z = s(0.7) + 1im*sign(s(0.9)-s(0.7))
	@test( closest(z,s) ≈ s(0.7) )
	@test( reflect(z,s) ≈ s(0.7) - (z-s(0.7)) )
	@test( all( abs(arg(s,s(t))-t) < 1e-11 for t in 0:0.1:1 ) )
end

@testset "Rays" begin
	s = Ray(Polar(2,0),pi/2)
	@test( isinf(arclength(s)) )
	@test( isleft(-1im,s) && !isleft(-1im,reverse(s)) )
	@test( real(s(0.23)) ≈ 2 )
	@test( imag(s(0.9)) > imag(s(0.7)) )
	@test( closest(5im,s) ≈ 2+5im )
	@test( all( abs(arg(s,s(t))-t) < 1e-11 for t in 0:0.1:1 ) )
	@test( angle(tangent(s,.1))≈π/2 )
	@test( tangent(s,0.1) ≈ CR.fdtangent(s,0.1) )
	@test( angle(tangent(reverse(s),1))≈-π/2  )
	s = Ray(Spherical(2im),pi,true)
	@test( imag(s(0.5)) ≈ 2 )
	@test( real(s(0.3)) < real(s(0.4)) )
	@test( tangent(s,0.6) ≈ CR.fdtangent(s,0.6) )
	@test( !isleft(4,s) && isleft(-1+3im,s) )
	@test( closest(-4+1im,s) ≈ -4+2im )
	@test( closest(6,s) ≈ 2im )
	@test( all( abs(arg(s,s(t))-t) < 1e-11 for t in 0:0.1:1 ) )
end

@testset "Intersections" begin 
	z = intersect(Circle(0,1),Circle(.2+0.5im,1.5))
	@test( all( @. abs(z-0)≈1 ) )
	@test( all( @. abs(z-(.2+0.5im))≈1.5 ) )
	z = intersect(Circle(0,1),Circle(1+1im,1.5))
	@test( all( @. abs(z-0)≈1 ) )
	@test( all( @. abs(z-(1+1im))≈1.5 ) )
	@test( isempty(intersect(Circle(0,1),Circle(.2+0.5im,.1))) )
	@test( isempty(intersect(Circle(0,1),Circle(.2+0.5im,6))) )

	z = intersect(Line(1,direction=1im),Line(-1,direction=1+1im))
	@test( z[1]≈1+2im )
	@test( isempty(intersect(Line(1,direction=1im),Line(-2,direction=1im))) )
	l = Line(2,direction=3+1im)
	@test( intersect(l,l+1e-15) ≈ l )

	z = intersect(Segment(0,1),Segment(.4-1im,.7+2im))
	@test( z[1]≈0.5 )
	z = intersect(Segment(2+3im,3+3im),Segment(-1+3im,4+3im))
	@test( z≈Segment(2+3im,3+3im) )
	z = intersect(2-Segment(2+3im,3+3im),2-Segment(3im,2.4+3im))
	@test( z≈2-Segment(2+3im,2.4+3im) )
	z = intersect(1im*Segment(2+3im,3+3im),1im*Segment(2.7+3im,3.4+3im))
	@test( z≈1im*Segment(2.7+3im,3+3im) )
	z = intersect(Segment(2+3im,3+3im),Segment(2.7+4im,3.4+6im))
	@test( isempty(z) )

	z = intersect(Ray(0,pi/4),Segment(.5-1im,.5+2im))
	@test( z[1]≈0.5+0.5im )
	@test( isempty(intersect(Ray(0,-3pi/4),Segment(.5-1im,.5+2im))) )
	z = intersect(Ray(0,pi/4),Ray(2+2im,pi/4))
	@test( z≈Ray(2+2im,pi/4) )
	z = intersect(Ray(0,pi/6),Ray(-1,-pi/2))
	@test( isempty(z) )

	z = intersect(3+Circle(0,1),3+Line(0.5,0.5+3im))
	@test( z[1]≈(3.5+sqrt(3)/2*1im) || z[1]≈(3.5-sqrt(3)/2*1im) )
	z = intersect(Circle(-2im,2),Line(3im,3))
	@test( isempty(z) )

	z = intersect(3+Ray(0.5+0.1im,pi/2,true),3+Circle(0,1))
	@test( z[1]≈(3.5+sqrt(3)/2*1im) && length(z)==1 )

	z = intersect(3+Circle(0,1),3+Segment(0.5,0.5+3im))
	@test( z[1]≈(3.5+sqrt(3)/2*1im) || z[1]≈(3.5-sqrt(3)/2*1im) )
	z = intersect(Circle(-2im,2),Segment(-1im,.5-1im))
	@test( isempty(z) )

end
	
@testset "Paths" begin
	S = Segment(1,1im)
	A = Arc(1im,-1+0.5im,-1)
	P = Path([S,A,-S])
	@test( all( point(P,[0,1,1.5,2.5,3]) .≈ [S(0),S(1),A(0.5),-S(0.5),-S(1)] ) )
	Q = 1 - 3im*P 
	@test( Q(1.5) ≈ 1 - 3im*A(0.5) )

	P = ClosedPath([S,1im*S,-S,-1im*S])
	@test( vertex(P,3) ≈ -1 )
	@test( arclength(P) ≈ 4*sqrt(2) )
	@test( isa(reverse(P),ClosedPath) )
	@test( all( point(P,[0,1,1.25,2.5,3,4]) .≈ [S(0),S(1),1im*S(.25),-S(0.5),-S(1),S(0)] ) )

end

@testset "Polygons" begin 
	s = Segment(2,2im) 
	p = Polygon([s,1im*s,-s,-1im*s]) 
	@test( arclength(p) ≈ 8*sqrt(2) ) 
	@test( angle(normal(p,1.1+length(p))) ≈ -π/4 )
	@test( winding(p,-0.4+0.5im) == 1 )
	@test( winding(reverse(p),-0.4+0.5im) == -1 )
	@test( winding(p,-4-0.5im) == 0 )
	@test( all( angles(p) .≈ 0.5*pi ) )
	@test( ispositive(p) )
	@test( !ispositive(reverse(p)) )

	p = Polygon([4,4+3im,3im,-2im,6-2im,6])
	@test( arclength(p) ≈ (3+4+5+6+2+2) )
	@test( tangent(p,2.3-length(p)) ≈ -5im )
	@test( winding(p,5-im) == 1 )
	@test( winding(p,-1) == 0 )
	@test( sum(angles(p)) ≈ 4*pi )

	p = CircularPolygon([Arc(1,2+1im,1im),Segment(1im,-1),Arc(-1,-0.5im,-1im),Segment(-1im,1)])
	@test( all(winding(p,z)==1 for z in [1+0.5im,1.7+1im,0,-1+0.05*exp(1im*pi/5),-1im+0.05*exp(1im*0.3*pi)]) )
	@test( all(winding(p,z)==0 for z in [-.999im,0.001-1im,-.999,-1.001,1.001,1.999im]) )
end

@testset "Unbounded polygons" begin
	p = Polygon([5,4+3im,3im,-2im,6-2im,(-pi/2,0)])
	a = angles(p)/pi
	@test( a[6] ≈ -0.5  )
	@test( sum(a.-1) ≈ -2 )

	p = Polygon([(pi/2,pi/2),5,4+3im,3im,-2im,6-2im])
	a = angles(p)/pi
	@test( abs(a[1]) < 1e-10 )
	@test( sum(a.-1) ≈ -2 )
	@test( all(winding(p,z)==1 for z in [1+2im,5-1im,5.5+6im]) )
	@test( all(winding(p,z)==0 for z in [-3im,3+5im,5.5-6im]) )
	
	p = Polygon([(-pi/2,pi/2),7,4+3im,3im,-2im,6-2im])
	a = angles(p)/pi
	@test( a[1] ≈ -1 )
	@test( sum(a.-1) ≈ -2 )
	@test( all(winding(p,z)==1 for z in [4,7-2im,9]) )
	@test( all(winding(p,z)==0 for z in [4+4im,6+2im,4-3im]) )
	
	p = Polygon([4+3im,7,(0,0),6-2im,-2im,3im])
	a = angles(p)/pi
	@test( a[3] ≈ -2 )
	@test( sum(a.-1) ≈ -2 )
end

@testset "Möbius" begin
	z = [0,1,2+2im]
	w = [Inf,-1im,-1]
	f = Möbius(z,w)
	@test( all(f.(z).≈w) )
	g = inv(f) 
	@test( all(g.(w).≈z) )
	c = Circle(z...)
	l = Line(w[2],w[3])
	@test( f(c)≈l && g(l)≈c )
	d = exterior(c) 
	h = halfplane(l) 
	@test( f(d)≈!h && g(h)≈!d )
	@test( f(!d)≈h && g(!h)≈d )
	f = Möbius([2,3+im,-5],z)
	g = Möbius(w,[2,3+im,-5])
	h = f ∘ g
	@test( all(h.(w).≈z) )
end