EXT_SUFFIX := $(python3-config --extension-suffix)

default: f90_functions$(EXT_SUFFIX)

f90_functions$(EXT_SUFFIX): f90_functions.f90 stats.o
	f2py3 -c -m f90_functions f90_functions.f90 stats.o
stats.o:stats.f90
	gfortran -c stats.f90