#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include "mesh.h"
#include "input.h"

void meshAlloc(grid * mesh){
    mesh->x = malloc(imax * sizeof(double));
}

void createMesh(grid * mesh){
    mesh->x[imax-1] = l;
    double dx = (double)((l/(imax-1))*2.0);
    for (int i = (imax-2); i >= 0; --i){
        mesh->x[i] = mesh->x[i+1] - dx;
    };
}

void freeMesh(grid * mesh){
    free(mesh->x);
}
