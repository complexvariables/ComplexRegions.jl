using ComplexValues,ComplexRegions
const Spherical = ComplexValues.Spherical
const Polar = ComplexValues.Polar

using Test
@testset "Utilities" begin
	 @test(ComplexRegions.@scaleto(1im,3im,[0.5,0.75]) ≈ [2.0im,2.5im])
	 @test(ComplexRegions.@scalefrom(1im,3im,[2im,1im,1.5im]) ≈ [0.5,0,0.25])
end

@testset "Binary functions on $S,$T" for S in [Polar,Spherical], T in [Polar,Spherical]
	s1 = Segment(3.0+5.0im,0)
	s2 = Segment(Spherical(1.0im),-1)
end 
	


