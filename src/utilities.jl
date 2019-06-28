macro scalefrom(a,b,t)
	return :(@. ($t-$a)/($b-$a))
end

macro scaleto(a,b,t)
	return :(@. $a + $t*($b-$a))
end

scalefrom(a,b,t) = @. (t-a)/(b-a)
scaleto(a,b,t) = @. a + t*(b-a)