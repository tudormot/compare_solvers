
CC=mpic++

#todo add flags
COMPFLAGS=-std=c++11 -Wall -Wextra -Weffc++

INCL_PETSC_HEADERS=-I/${PETSC_DIR}/include -I/${PETSC_DIR}/${PETSC_ARCH}/include

#note:make sure PETSC_DIR env variable exists. If not, this makefile fails
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

#this rule makes it so that everything is recompiled when a header file changes(provided DEPS(dependencies) var contains all headers)
DEPS=lin_sys.h file_reader.h
%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(COMPFLAGS)

debug: main.cpp lin_sys.cpp file_reader.cpp
	${CC} -o debug ${COMPFLAGS} main.cpp ${INCL_PETSC_HEADERS} ${PETSC_LIB}

#TODO:create rule for release build

