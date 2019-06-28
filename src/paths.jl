abstract type AbstractPath end

struct Path <: AbstractPath 
	curve
	offsets
	function Path(p::AbstractVector;tol::Real=1e-13)
		for k = 1:length(p)-1
			@assert p[k] isa AbstractCurve
			@assert isapprox(point(p[k],1.0),point(p[k+1],0.0),rtol=tol,atol=tol) "Curve endpoints do not match for pieces $(k) and $(k+1)"
		end
		new(p,cumsum([arclength(c) for c in p]))
	end
end
Path(c::AbstractCurve) = Path([c])

arclength(p::Path) = p.offsets[end]

function point(p::Path,t::Real)
	#@assert (0 ≤ t ≤ 1) "Parameter is out of the range [0,1]."
	offset = p.offsets/p.offsets[end]
	j = findlast(t .> offset)
	if isnothing(j)
		s = scalefrom(0,offset[1],t) 
		p.curve[1](s)
	else
		s = scalefrom(offset[j],offset[j+1],t) 
		p.curve[j+1](s)
	end
end

(p::Path)(t::Real) = point(p,t)

iterate(p::Path,state=1) = state > length(p.curve) ? nothing : (p[state], state+1)
eltype(::Type{Path}) = AbstractCurve 
length(p::Path) = length(p.curve)
getindex(p::Path,k) = getindex(p.curve,k)