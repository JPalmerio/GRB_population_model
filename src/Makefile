EXT_SUFFIX := $(python3-config --extension-suffix)

default: f90_functions$(EXT_SUFFIX)

f90_functions$(EXT_SUFFIX): f90_functions.f90 stats_f90.o
	@echo 'Extension suffix is: $(EXT_SUFFIX)'
	f2py3 -c -m f90_functions f90_functions.f90 stats_f90.o
stats_f90.o:stats_f90.f90
	gfortran -c stats_f90.f90