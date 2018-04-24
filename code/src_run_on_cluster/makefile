
COMPFLAGS= -std=c++11

#all flags that have to do with pardiso library

#note:make sure PETSC_DIR env variable exists. If not, this makefile fails
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
INCL_PETSC_HEADERS=-I/${PETSC_DIR}/include -I/${PETSC_DIR}/${PETSC_ARCH}/include
#CXX= set by PETSC include

#this rule makes it so that everything is recompiled when a header file changes(provided DEPS(dependencies) var contains all headers)
DEPS=pardiso_solver.h petsc_solver.h sys_mat.h sys_solver.h testing.h timer.h
OBJ=main.o intel_pardiso_solver.o petsc_solver.o sys_mat.o testing.o timer.o 

%.o: %.cpp $(DEPS)
	$(PETSC_CXXCOMPILE_SINGLE) $< $(COMPFLAGS) ${INCL_PETSC_HEADERS} ${MKL_INC}
#$(PETSC_CXXCOMPILE) -o $@ $< $(COMPFLAGS)

release: $(OBJ)
	${CLINKER} -o release $(OBJ) ${MKL_SHLIB} ${PETSC_LIB}



