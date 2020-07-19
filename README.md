# mathmodules

A collection of modules for building mathematical software (in particular numerical physics programs), mostly in fortran.

## Make target: quadrature.o

Numerical integration module (fortran).
Provides routines for computing Gaussian quadrature nodes and weights
(refers to 'normal' Gaussian quadrature as Gauss-Legendre quadrature), 
and Gauss-Legendre-Lobatto nodes and weights. 
Also provides routines for computing Legendre polynomials and their derivatives,
as well as Lagrange interpolating polynomials. 
All the implementations are fairly "naive" in the sense that they're just coded straight from
the mathematical definitions. This is fine for my use since I typically call these methods only at the start
of a program, and store the results I need, but if you need very fast implementations, perhaps look elsewhere.
Read the comments in the source file for more information.

## Make target: pymodule_quadrature

Python wrapper for the above (uses f2py3, so you must have that installed, e.g. via pip3).
Once buit, you can import the module from python and call its functions.
For more information, look up f2py.

## Make target: geometry_quad4.o

A collection of useful geometry functions for dealing with coordinate mappings from cartesian to deformed quadrilaterals (fortran).
Useful for 2D finite element solvers and similar.
I use this for a 2D quadrilateral mesh system which is not yet included in this repository. 
Read the comments in the source file for more information.

## Make target: pymodule_geometry_quad4

Python wrapper for the above (uses f2py3, so you must have that installed, e.g. via pip3).
Once buit, you can import the module from python and call its functions.
For more information, look up f2py.

