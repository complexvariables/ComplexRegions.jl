module ComplexRegionsPythonCallExt

using ComplexRegions, PythonCall

function ComplexRegions.Curve(f::Py, dfdt::Py,  args...; kw...)
    tangent = t -> pyconvert(Number, pycall(dfdt, t))
    ComplexRegions.Curve(t -> pyconvert(Number, pycall(f, t)), tangent, args...; kw...)
end

function ComplexRegions.Curve(f::Py, args...; kw...)
    ComplexRegions.Curve(t -> pyconvert(Number, pycall(f, t)), args...; kw...)
end

function ComplexRegions.ClosedCurve(f::Py, dfdt::Py,  args...; kw...)
    tangent = t -> pyconvert(Number, pycall(dfdt, t))
    ComplexRegions.ClosedCurve(t -> pyconvert(Number, pycall(f, t)), tangent, args...; kw...)
end

function ComplexRegions.ClosedCurve(f::Py, args...; kw...)
    ComplexRegions.ClosedCurve(t -> pyconvert(Number, pycall(f, t)), args...; kw...)
end

end  # module
