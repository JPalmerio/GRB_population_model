EXT_SUFFIX := $(python3-config --extension-suffix)

default: f90_functions$(EXT_SUFFIX)


stats_f90.o:stats_f90.f90
	gfortran -c stats_f90.f90

f90_functions$(EXT_SUFFIX): f90_functions.f90 stats_f90.o stats_f90.mod
	f2py3 -c -m f90_functions f90_functions.f90 stats_f90.o
clean:
	rm -rf *.mod *.o
	rm -rf f90_functions.cpython-37m-darwin.so
	rm -rf f90_functions.cpython-37m-darwin.so.dSYM