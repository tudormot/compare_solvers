
#CC=mpic++ set by petsc 

COMPFLAGS=-std=c++11 -Wall -Wextra -Weffc++

#all flags that have to do with pardiso library
#TODO:currently assumed pardiso library placed in /usr/lib..Should change that to make it independent
PARDISO_LIB=-lpardiso500-GNU481-X86-64 -llapack  -lf77blas -lcblas -latlas -lgfortran -fopenmp -lpthread -lm

PETSC_DIR=/home/mot_tudor/workspace/unijob/various_libs/petsc-3.9.0
PETSC_ARCH=arch-linux2-c-debug
INCL_PETSC_HEADERS=-I/${PETSC_DIR}/include -I/${PETSC_DIR}/${PETSC_ARCH}/include

#note:make sure PETSC_DIR env variable exists. If not, this makefile fails
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

#this rule makes it so that everything is recompiled when a header file changes(provided DEPS(dependencies) var contains all headers)
DEPS=pardiso_solver.h petsc_solver.h sys_mat.h sys_solver.h testing.h timer.h
%.o: %.cpp $(DEPS)
	$(CXX) -c -o $@ $< $(COMPFLAGS)

Debug: main.cpp pardiso_solver.cpp petsc_solver.cpp sys_mat.cpp testing.cpp timer.cpp
	${CXX} -o debug ${COMPFLAGS} main.cpp pardiso_solver.cpp petsc_solver.cpp sys_mat.cpp testing.cpp timer.cpp ${INCL_PETSC_HEADERS} ${PARDISO_LIB} ${PETSC_LIB}

#TODO:create rule for release build

