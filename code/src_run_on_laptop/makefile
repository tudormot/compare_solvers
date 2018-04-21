
#CC=mpic++ set by petsc 

COMPFLAGS=-std=c++11 -Wall -Wextra 

#all flags that have to do with pardiso library
#TODO:currently assumed pardiso library placed in /usr/lib..Should change that to make it independent
PARDISO_LIB=-lpardiso500-GNU481-X86-64 -llapack  -lf77blas -lcblas -latlas -lgfortran -fopenmp -lpthread -lm

PETSC_DIR=/home/mot_tudor/workspace/unijob/various_libs/petsc-3.9.0
PETSC_ARCH=arch-linux2-c-debug
INCL_PETSC_HEADERS=-I/${PETSC_DIR}/include -I/${PETSC_DIR}/${PETSC_ARCH}/include

#note:make sure PETSC_DIR env variable exists. If not, this makefile fails
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

#CXX= set by PETSC include

#this rule makes it so that everything is recompiled when a header file changes(provided DEPS(dependencies) var contains all headers)
DEPS=pardiso_solver.h petsc_solver.h sys_mat.h sys_solver.h testing.h timer.h
OBJ=main.o pardiso_solver.o petsc_solver.o sys_mat.o testing.o timer.o
#ODIR = obj #used for putting object files in a separate directory
%.o: %.cpp $(DEPS)
	$(PETSC_CXXCOMPILE_SINGLE) $< $(COMPFLAGS) ${INCL_PETSC_HEADERS}

release: $(OBJ)
	${CLINKER} $(OBJ) ${PARDISO_LIB} ${PETSC_LIB}
