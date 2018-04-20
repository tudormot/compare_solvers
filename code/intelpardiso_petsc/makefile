
COMPFLAGS=-std=c++11 -O3

#all flags that have to do with pardiso library
INTELPARDISO_LIB=-lmkl_rt

#note:make sure PETSC_DIR env variable exists. If not, this makefile fails
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
INCL_PETSC_HEADERS=-I/${PETSC_DIR}/include -I/${PETSC_DIR}/${PETSC_ARCH}/include
#CC= set by PETSC include

#this rule makes it so that everything is recompiled when a header file changes(provided DEPS(dependencies) var contains all headers)
DEPS=pardiso_solver.h petsc_solver.h sys_mat.h sys_solver.h testing.h timer.h
%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(COMPFLAGS)

release: main.cpp pardiso_solver.cpp petsc_solver.cpp sys_mat.cpp testing.cpp timer.cpp
	${CC} -o time_solvers ${COMPFLAGS} main.cpp pardiso_solver.cpp petsc_solver.cpp sys_mat.cpp testing.cpp timer.cpp ${INCL_PETSC_HEADERS} ${INTELPARDISO_LIB} ${PETSC_LIB}



