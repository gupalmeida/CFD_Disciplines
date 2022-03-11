# Makefile for compiling C++ code

CC=gcc

CFLAGS=-O

PROJ_NAME=shocktube

C_SOURCE=$(wildcard *.c)
H_SOURCE=$(wildcard *.h)

# object files
#OBJ=$(subst .C,.o,$(subst source,objects,$(C_SOURCE)))

$(PROJ_NAME): $(C_SOURCE) $(H_SOURCE)
	$(CC) $(CFLAGS) $(C_SOURCE) $(H_SOURCE) -lm -o $(PROJ_NAME) 

clean:
	@ echo "Cleaning object files\n"
	@ rm -rf *.o 

