CFLAGS=-c -Wall -pedantic
CFLAGS+=-O3 -march=native -funroll-loops -ffast-math -Wall -Wno-tabs -pedantic
#CFLAGS+=-g

all: quadrature.o pymodule_quadrature geometry_quad4.o pymodule_geometry_quad4

quadrature.o: quadrature.f90
	gfortran $(CFLAGS) $^ -o $@ 

pymodule_quadrature: quadrature.f90
	f2py3 -c $^ -m quadrature 

# You'll have to link LAPACK when you build something with this module
geometry_quad4.o: geometry_quad4.f90
	gfortran $(CFLAGS) $^ -o $@ 

pymodule_geometry_quad4: geometry_quad4.f90
	f2py3 -c $^ -m geometry_quad4 --link-lapack


clean:
	rm -f *.so *.o *.mod

.PHONY: pymodule_quadrature pymodule_geometry_quad4 all clean
