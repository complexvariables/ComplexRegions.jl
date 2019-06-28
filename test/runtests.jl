using ComplexValues,ComplexRegions

using Test
@testset "Utilities" begin
	 @test(ComplexRegions.@scaleto(1im,3im,[0.5,0.75]) ≈ [2.0im,2.5im])
	 @test(ComplexRegions.@scalefrom(1im,3im,[2im,1im,1.5im]) ≈ [0.5,0,0.25])
end

@testset "Curves" begin
	c = Circle(Spherical(1-1im),sqrt(2))
	@test( arclength(c) ≈ 2*sqrt(2)*pi )
	c = Circle(1-1im,sqrt(2))
	@test( point(c,.25) ≈ complex(1,sqrt(2)-1) )
	c = Circle(1f0,-1im,0)
	@test( c.radius ≈ 1/sqrt(2f0) )
	zz = point(5im - c/3im,.23)
	@test( abs(zz-(5im-c.center/3im)) ≈ c.radius/3 )

	@test( Circle(1+3im,Polar(1-1im),1.0) isa Line )

	@test( Line(1,5) isa Line )
	l = Line(1im,direction=1+2im)
	dz = point(l,.6)-point(l,.1)
	@test( angle(dz) ≈ angle(1+2im) )
	zz = point(5im - l/3im,.23)
	z0 = 5im - l.base/3im
	@test( angle(zz-z0) ≈ angle(l.direction/3im) )

	a = Arc(1.0,1im,center=0)
	zz = 1/sqrt(2)*(1+1im)
	@test( point(a,0.5) ≈ zz)
	b = Arc(point(a,0),point(a,.25),point(a,1))
	@test( point(b,0.5) ≈ zz)
	@test( point(0.1 - 3im*a,0.5) ≈ 0.1 - 3im*zz)

	#s1 = Segment(3.0+5.0im,0)
	#s2 = Segment(Spherical(1.0im),-1)

end 
	


