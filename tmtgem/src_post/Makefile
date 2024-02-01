# Makefile for the use of gfotran T. MINAMI on 2020.12.10
# slice.exe is added on 2017.12.20
# Coded on August 15, 2015 by T. MINAMI

FC=gfortran
MPIFC=mpif90
#OPT=-i-dynamic
VPATH = ./ ../common
#OPT0= -O3 -xAVX -mkl=parallel
#OPT0= -O -fpe0 -check bounds -traceback -g -mkl=parallel

#2016.03.03 solver with many modules

target2=slice.exe
target6=change_model2cond.exe

all: $(target2) $(target6)

objs2= m_constants.o  m_outerinnerproduct.o m_matrix_gfort.o m_fem_util.o m_mesh_type.o  sort_index.o m_param.o m_post_param.o m_slicesub.o slice.o

objs6=m_param.o m_fem_util.o m_outerinnerproduct.o sort_index.o m_cond_type.o m_mesh_type.o m_triangle.o m_cond_model_type.o change_model2cond.o



$(target2) : $(objs2)
	$(FC) $(OPT0) $(objs2) -o $(target2)

$(target6) : $(objs6)
	$(FC) $(OPT0) $(objs6) -o $(target6)

# rules
.SUFFIXES : .f90 .o
.f90.o:
	$(FC) -c $(OPT0) $< -o $@
.f.o:
	$(FC) -c $(OPT0) $< -o $@

clean:
	rm -f *.o    
	rm -f *.mod  
	rm -f *.exe  
