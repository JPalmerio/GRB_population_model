# Compilateur

#F90=gfortran -ffree-line-length-0
F90=mpif90 -fcheck='bounds' -ffree-line-length-0 -g -fbacktrace -ffpe-trap=zero,overflow,underflow

GRB_population_MC.exe:mtmod.o stats.o physics.o cosmology.o param.o GRB_population_MC.o
	$(F90) -o GRB_population_MC.exe mtmod.o stats.o physics.o cosmology.o param.o GRB_population_MC.o
mtmod.o:mtmod.f90
	$(F90) -c mtmod.f90
stats.o:stats.f90
	$(F90) -c stats.f90
physics.o:physics.f90
	$(F90) -c physics.f90
cosmology.o:cosmology.f90
	$(F90) -c cosmology.f90
param.o:param.f90
	$(F90) -c param.f90
GRB_population_MC.o:GRB_population_MC.f90
	$(F90) -c GRB_population_MC.f90

clean:
	rm *.o ; rm *.mod