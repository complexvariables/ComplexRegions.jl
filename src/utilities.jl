# Scale from/to [0,1] 
scalefrom(a,b,t) = @. (t-a)/(b-a)
scaleto(a,b,t) = @. a + t*(b-a)

#inf(α) = Polar(Inf,α*π)

# Unique real roots of a quadratic.
function realroots(a,b,c)
	a==0 && return [-c/b]
	bh = -b/2
	d = bh^2 - a*c 
	if d < 0
		[]
	elseif d==0
		[bh/a]
	else
		r = (bh + sign(bh)*sqrt(d)) / a 
		[r,c/(r*a)]
	end
end  

# Are three points arranged in counterclockwise order? 
function isccw(a::Number,b::Number,c::Number)
	v(z) = SVector(1,real(z),imag(z))
	det([v(a) v(b) v(c)]) > 0
end
