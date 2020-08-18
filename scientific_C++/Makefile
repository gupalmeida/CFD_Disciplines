# Makefile for compiling C++ code
CC=g++

CFLAGS=-g

PROJ_NAME=fdSolver

#SRCS=$(wildcard ./source/*.C)
#H_SOURCE=$(wildcard ./source/*.H)

# object files
#OBJ=$(subst .C,.o,$(subst source,objects,$(C_SOURCE)))
all: fdSolver clean

fdSolver: main.o IOTools.o Vector.o Matrix.o LinSys.o FDGrid.o
	$(CC) $(CFLAGS) main.o IOTools.o Vector.o Matrix.o LinSys.o FDGrid.o BoundaryConditions.o BvpODE.o -o $(PROJ_NAME) 

IOTools.o: IOTools.C IOTools.H tmpIOTools.H
	$(CC) $(CFLAGS) -c IOTools.C

Vector.o: Vector.C Vector.H
	$(CC) $(CFLAGS) -c Vector.C

Matrix.o: Matrix.C Matrix.H Vector.o
	$(CC) $(CFLAGS) -c Matrix.C

LinAlg.o: LinSys.C LinSys.H Matrix.o Vector.o
	$(CC) $(CFLAGS) -c LinSys.C

FDGrid.o: FDGrid.C FDGrid.H Node.H
	$(CC) $(CFLAGS) -c FDGrid.C

BondaryConditions.o: BoundaryConditions.C BoundaryConditions.H
	$(CC) $(CFLAGS) -c BoundaryConditions.C

BvpODE.o: BvpODE.C BvpODE.H BoundaryConditions.o ScndOrderODE.H FDGrid.o
	$(CC) $(CFLAGS) -c BvpODE.C

main.o: main.C IOTools.o BvpODE.o
	$(CC) $(CFLAGS) -c main.C

clean:
	@ clear
	@ echo "\n======================\n"
	@ echo "Cleaning object files\n"
	@ echo "======================\n"
	@ rm -rf *.o 

