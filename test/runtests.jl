using ComplexValues,ComplexRegions

using Test
@testset "Utilities" begin
	 @test(ComplexRegions.@scaleto(1im,3im,[0.5,0.75]) ≈ [2.0im,2.5im])
	 @test(ComplexRegions.@scalefrom(1im,3im,[2im,1im,1.5im]) ≈ [0.5,0,0.25])
end

@testset "Curves" begin
	@test( Curve(t -> 2*cos(t) + 1im*sin(t),-1,1) isa Curve)
	@test( ClosedCurve(t -> 2*cos(t) + 1im*sin(t),0,2π) isa ClosedCurve)
end

@testset "CirclesArcs" begin
	c = Circle(Spherical(1-1im),sqrt(2))
	@test( arclength(c) ≈ 2*sqrt(2)*pi )
	@test( dist(-1+1im,c) ≈ sqrt(2) )
	@test( closest(1+4im,c) ≈ 1+1im*(sqrt(2)-1) )
	c = Circle(1-1im,sqrt(2))
	@test( point(c,.25) ≈ complex(1,sqrt(2)-1) )
	@test( 2/c isa Line )
	c = Circle(1f0,-1im,0)
	@test( c.radius ≈ 1/sqrt(2f0) )
	zz = point(5im - c/3im,.23)
	@test( abs(zz-(5im-c.center/3im)) ≈ c.radius/3 )
	@test( Circle(1+3im,Polar(1-1im),1.0) isa Line )

	a = Arc(1.0,1im,center=0)
	zz = 1/sqrt(2)*(1+1im)
	@test( point(a,0.5) ≈ zz)
	b = Arc(point(a,0),point(a,.25),point(a,1))
	@test( point(b,0.5) ≈ zz)
	@test( point(0.1 - 3im*a,0.5) ≈ 0.1 - 3im*zz)
	@test( closest(-1+5im,a) ≈ 1im )
	@test( closest(2-5im,a+1) ≈ 1+1 )
	@test( dist(3im+.5*2im*exp(1im*pi/5),3im+2im*a) ≈ 2*0.5 )
end 

@testset "LinesSegments" begin
	@test( Line(1,5) isa Line )
	l = Line(1im,direction=1+2im)
	dz = point(l,.6)-point(l,.1)
	@test( angle(dz) ≈ angle(1+2im) )
	@test( 1/l isa Circle )
	zz = point(5im - l/3im,.23)
	z0 = 5im - l.base/3im
	@test( angle(zz-z0) ≈ angle(l.direction/3im) )
	z = l(0.3) + 1im*sign(l.direction)
	@test( dist(z,l) ≈ 1 )
	@test( closest(z,l) ≈ l(0.3) )

	s = Segment(1,3.0+5.0im)
	zz = 2 + 2.5im
	@test( point(2 - 3im*s,0.5) ≈ 2 - 3im*zz)
	@test( closest(4+6im,s) ≈ 3+5im )
	@test( dist(-1,s) ≈ 2 )
	z = s(0.7) + 1im*sign(s.line.direction)
	@test( closest(z,s) ≈ s(0.7) )
end
	
@testset "Paths" begin
	S = Segment(1,1im)
	A = Arc(1im,-1+0.5im,-1)
	l1,l2 = arclength(S),arclength(A)
	l = l1+l2
	P = Path([S,A,-S])
	t = (l1+0.5l2)/(2l1+l2)
	@test( P(t) ≈ A(0.5) )
	Q = 1 - 3im*P 
	@test( Q(t) ≈ 1 - 3im*A(0.5) )

	P = ClosedPath([S,1im*S,-S,-1im*S])
	@test( vertex(P,3) ≈ -1 )
	@test( arclength(P) ≈ 4*sqrt(2) )
end

