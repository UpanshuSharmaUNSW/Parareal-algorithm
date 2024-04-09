# this variable defines the compiler's name
CC=g++

# path to libraries (external)
LDFLAGS=-L/usr/lib/
# LIBRARY=-lstdc++  -lm -lblas -llapack  -lg2c -lgslcblas

# place where
INCLUDE_PATH=-I./

# optimizations
CPPFLAGS= -Wall -O2 -ffast-math -funroll-loops ${INCLUDE_PATH}

# conditional compilation
# DEFINES = -DHAVE_LAPACK -DHAVE_CBLAS
DEFINES=

# rule for compiling each of the cpp files
.cpp.o:
	${CC} ${CPPFLAGS} ${DEFINES} ${INCLUDE_PATH} -g -c $*.cpp

OBJ= gaussien.o input.o hamiltonian.o matrice_GS.o algorithm_natoms.o vecteur_FL.o particle.o main.o

# takes all the *.o files and combines them as the CG++ executable
all: ${OBJ}
	${CC} ${CPPFLAGS} ${DEFINES}  ${OBJ} -o CG++ ${LDFLAGS}

# cleaning
clean:
	/bin/rm -rf *~ *.o CG++
