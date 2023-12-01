##############################################################
# Makefile: Cornell Multigrid Coupled Tsunami (COMCOT) Model #
# LAST REVISED: XIAOMING WANG  30/10/2008                    #
##############################################################

#-------- FOR WINDOWS ---------------------------
#F90         =  ifort
#OBJ         =  comcot.exe
#SRC         =  comcot.f90 type_module.f90 initialization.f90 output.f90 deform.f90 mass.f90 moment.f90 boundaries.f90 all_grids.f90 hotstart.f90 landslide.f90 wavemaker.f90 dispersion.f90
#OBJECTS     =  comcot.obj type_module.obj initialization.obj output.obj deform.obj mass.obj moment.obj boundaries.obj all_grids.obj hotstart.obj landslide.obj wavemaker.obj dispersion.obj
#MODULES     =  fault_params.mod landslide_params.mod layer_params.mod wave_params.mod bci_params.mod

#OPT          = /O2 /F1000000000
#OPT         = -mkl=parallel -g -traceback

#$(OBJ):  $(SRC)
#	$(F90) $(SRC) $(OPT) -o $(OBJ)
	
#CLEAN:
#	RM $(OBJECTS)
#	RM $(MODULES)
	
##--------- FOR LINUX ---------------------

F90         =  ifort
OBJ         =  comcot.exe
OBJECTS     =  type_module.o initialization.o output.o deform.o mass.o moment.o boundaries.o all_grids.o hotstart.o landslide.o wavemaker.o dispersion.o comcot.o

#OPT         =  -O2 -F1000000000
OPT = -qmkl=parallel -O -g -traceback

$(OBJ):  $(OBJECTS)
	$(F90) $(OBJECTS) $(OPT) -o $(OBJ)

.SUFFIXES : .f90 .o
.f90.o:
	$(F90) -c $(OPT) $< -o $@
.f.o:
	$(F90) -c $< -o $@
#
clean:
	rm *.o
	rm *.mod
	rm *.exe

