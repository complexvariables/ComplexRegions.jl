macro scalefrom(a,b,t)
	return :(@. ($t-$a)/($b-$a))
end

macro scaleto(a,b,t)
	return :(@. $a + $t*($b-$a))
end

scalefrom(a,b,t) = @. (t-a)/(b-a)
scaleto(a,b,t) = @. a + t*(b-a)

inf(α) = Polar(Inf,α*π)

function quadroots(a,b,c)
	a==0 && return [-c/b]
	bh = -b/2
	d = bh^2 - a*c 
	if d==0
		return [bh/a,bh/a]
	else
		r = (bh + sign(bh)*sqrt(d))/a 
		return [r,c/(r*a)]
	end
end  

function isccw(a::Number,b::Number,c::Number)
	v(z) = SVector(1,real(z),imag(z))
	det([v(a) v(b) v(c)]) > 0
end
