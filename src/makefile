#this is an incomplete makefile, that has to be modified in order to compile the program.
#As specified in the README.md file, three variables have to be set correctly in order
#to succesfully compile the program: MKL_INC, MKL_SHLIB, PETSC_DIR, PETSC_ARCH

#these should be completed, depending on your MKL and PETSc installations:
MKL_INC=
MKL_SHLIB=
PETSC_DIR=
PETSC_ARCH= 


#the following lines present an example on how to set these variables(based on my laptop's settings):
#INTELROOT = /opt/intel
#MKLROOT= /opt/intel/mkl
#MKL_INC=  -m64 -I${MKLROOT}/include
#MKL_SHLIB= -L${MKLROOT}/lib/intel64 -L${INTELROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl
#PETSC_DIR=/home/mot_tudor/workspace/unijob/various_libs/petsc-3.9.0
#PETSC_ARCH=arch-linux2-c-debug

#note:make sure PETSC_DIR env variable exists. If not, this makefile fails
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

COMPFLAGS=-std=c++11 -Wall -Wextra

DEPS=intel_pardiso_solver.h petsc_solver.h sys_mat.h sys_solver.h testing.h timer.h
OBJ=main.o intel_pardiso_solver.o petsc_solver.o sys_mat.o testing.o timer.o 

%.o: %.cpp $(DEPS)
	$(PETSC_CXXCOMPILE_SINGLE) $< $(COMPFLAGS) ${INCL_PETSC_HEADERS} ${MKL_INC}
	
all: $(OBJ)
	${CLINKER} -o time_solvers  $(OBJ) ${MKL_SHLIB} ${PETSC_LIB}


