#ifndef MESH_H
#define MESH_H

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include "input.h"

typedef struct grid{
    /* x is a unidimensional array
    * used to allocate x-coordinates
    * of the computational mesh */
    
    double *x;
}grid;

void meshAlloc(grid *);
void createMesh(grid *);
void freeMesh(grid *);

#endif
