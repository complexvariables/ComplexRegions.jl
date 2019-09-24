# Documentation for methods defined multiple times in parallel ways, to cut down on repetition/duplication when requesting help.

@doc """
	Complex(::AbstractCurve) 
Interpret a curve as having points of type Complex.
""" Complex(::AbstractCurve) 
@doc """
	Polar(::AbstractCurve) 
Interpret a curve as having points of type Polar.
""" Polar(::AbstractCurve) 
@doc """
	Spherical(::AbstractCurve) 
Interpret a curve as having points of type Spherical.
""" Spherical(::AbstractCurve) 

@doc """
	winding(P,z)
Compute the winding number of a closed curve or path `P` about the point `z`. Each counterclockwise rotation about `z` contributes +1, and each clockwise rotation about it counts -1. The winding number is zero for points not in the region enclosed by `P`. 

The result is unreliable for points lying on `P`, for which the problem is ill-posed.
""" winding(::Union{AbstractClosedCurve,AbstractClosedPath})

# curve/path/region
AbstractCPR = Union{AbstractCurve,AbstractPath,AbstractRegion}
@doc """
	X + z
	z + X 
Translate a curve, path, or region `X` by a complex number `z`. 
""" +(::AbstractCPR)

@doc """
	-X 
Negate a curve, path, or region `X` (reflect through the origin).
""" -(::AbstractCPR)

@doc """
	X - z
Translate a curve, path, or region `X` by a number `-z`.
""" -(::AbstractCPR,::Number)

@doc """
	z - X 
Negate a curve, path, or region `X` (reflect through the origin) and translate by `z`.
""" -(::Number,::AbstractCPR)

@doc """
	z*X 
	X*z 
Multiply the curve, path, or region `X` by complex number `z`; i.e., scale and rotate it about the origin.
""" *(::AbstractCPR)

@doc """
	X/z 
Multiply the curve, path, or region `X` by the number `1/z`; i.e., scale and rotate it about the origin.
""" /(::AbstractCPR,::Number)

@doc """
	z/X 
Invert the curve, path, or region `X` pointwise and multiply by the number `z`.
""" /(::Number,::AbstractCPR)

@doc """
	inv(X)
Invert a curve, path, or region pointwise.	
""" inv(::AbstractCPR)

@doc """
	conj(X) 
Construct the complex conjugate of curve, path, or region `X`. (Reverses the orientation of a curve or path.)
""" conj(::AbstractCPR)

# curve/path 
AbstractCP = Union{AbstractCurve,AbstractPath}

@doc """ 
	reverse(X) 
Construct a curve or path identical to `X` except with opposite direction of parameterization.
""" reverse(::AbstractCP)

@doc """ 
	arclength(X) 
Fetch or compute the arc length of the curve or path `X`.

# Example
```
julia> ellipse = ClosedCurve( t->cos(t)+2im*sin(t), 0,2π );
julia> arclength(ellipse)  # good to about 10 digits
9.688448219981513
```
""" arclength(::AbstractCP)

AbstractCCP = Union{AbstractClosedCurve,AbstractClosedPath}

@doc """
	isinside(X,z)
Detect whether `z` lies inside the closed curve or path `X`. For a bounded path, this always means the bounded region enclosed by the curve, regardless of orientation; for an unbounded path, it means the region "to the left" as one walks along the path.
""" isinside(z::Number,C::AbstractCCP)

@doc """
	isoutside(X,z)
Detect whether `z` lies outside the closed curve or path `X`. For a bounded path, this always means the unbounded region complementary to the enclosure of the curve, regardless of orientation; for an unbounded path, it means the region "to the right" as one walks along the path.
""" isoutside(z::Number,C::AbstractCCP)