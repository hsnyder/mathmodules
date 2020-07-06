# mathmodules

A collection of modules for building mathematical software, mostly in fortran.

## Make target: quadrature.o

Numerical integration module (fortran).
Read the comments in the source file for more information.

## Make target: pymodule_quadrature

Python wrapper for the above (uses f2py3, so you must have that installed, e.g. via pip3).
Once buit, you can import the module from python and call its functions.
For more information, look up f2py.

## Make target: geometry_quad4.o

A collection of useful geometry functions for dealing with coordinate mappings from cartesian to deformed quadrilaterals (fortran).
Useful for 2D finite element solvers and similar.
Read the comments in the source file for more information.

## Make target: pymodule_geometry_quad4

Python wrapper for the above (uses f2py3, so you must have that installed, e.g. via pip3).
Once buit, you can import the module from python and call its functions.
For more information, look up f2py.

