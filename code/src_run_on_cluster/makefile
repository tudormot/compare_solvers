
COMPFLAGS= -std=c++11 -g

#all flags that have to do with pardiso library. On cluster these are not required as the environmental variables take care of it
laptop: INTELROOT = /opt/intel
laptop: MKLROOT= /opt/intel/mkl
laptop: MKL_INC=  -m64 -I${MKLROOT}/include
laptop: MKL_SHLIB= -L${MKLROOT}/lib/intel64 -L${INTELROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl

#note:make sure PETSC_DIR env variable exists. If not, this makefile fails
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
INCL_PETSC_HEADERS=-I/${PETSC_DIR}/include -I/${PETSC_DIR}/${PETSC_ARCH}/include

#this rule makes it so that everything is recompiled when a header file changes(provided DEPS(dependencies) var contains all headers)
DEPS=pardiso_solver.h petsc_solver.h sys_mat.h sys_solver.h testing.h timer.h
OBJ=main.o intel_pardiso_solver.o petsc_solver.o sys_mat.o testing.o timer.o 

%.o: %.cpp $(DEPS)
	$(PETSC_CXXCOMPILE_SINGLE) $< $(COMPFLAGS) ${INCL_PETSC_HEADERS} ${MKL_INC}
#$(PETSC_CXXCOMPILE) -o $@ $< $(COMPFLAGS)

laptop: $(OBJ)
	${CLINKER} -o laptop  $(OBJ) ${MKL_SHLIB} ${PETSC_LIB}
release: $(OBJ)
	${CLINKER} -o release $(OBJ) ${MKL_SHLIB} ${PETSC_LIB}



