# Usage from Python

You can call the functions in this package from Python using the [`PythonCall`/`JuliaCall`](https://juliapy.github.io/PythonCall.jl/stable/) package. 

## Installation

It's recommended to create a new virtual environment to try this out. In Python, you need to install `juliacall` via

```bash
pip install juliacall
```

Then, start Python and run:

```python
from juliacall import Main as jl
```

This will download and initialize a copy of Julia. Finally, you need to install this package in that Julia environment:

```python
jl.seval('using Pkg; Pkg.add("ComplexRegions")')
```

That should be all you need to set up in the Python environment.

## Usage

In each new Python session, you need to load the packages as follows:

```python
from juliacall import Main as jl
jl.seval('using ComplexRegions, PythonCall')
```

All the functions and constants exposed to Julia by this package are available using the `jl` object. For example, to use the discrete AAA algorithm:

```python
import numpy as np    # if installed in Python
a = jl.Arc(-1, 1, -1j)
print(a)
```

```
Julia:
Arc{Float64} in the complex plane:
   fraction 0.75 of (Circle(0.0+0.0im,1.0,cw)) starting at 0.5
```

Above, `a` is a wrapped Julia object that you can use in Python. For example, you can find a point on the arc:

```python
a(0.5)
```

```
(0.7071067811865476+0.7071067811865475j)
```

You can get information about the approximation using any documented function in the package, e.g.:

```python
print(jl.inv(a)) 
```

```
Julia:
Arc{Float64} in the complex plane:
   fraction 0.75 of (Circle(0.0+0.0im,1.0,ccw)) starting at 0.5
```

``` python
print(jl.arclength(a))
```

```
4.71238898038469
```

## Native Python interface

There is an alternative (currently experimental) native Python interface called `cxregions`. See the [GitHub home](https://github.com/complexvariables/cxregions) for details.
