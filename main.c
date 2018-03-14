#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "input.h"
#include "mesh.h"
#include "aux.h"

void writeSolution(results * , grid * );

int main(){
    
    /* parameter definition */
    grid mesh;
    results solution;
    
    /* generating mesh */
    meshAlloc(&mesh);
    createMesh(&mesh);
    
    /* initializing variables */
    allocSolution(&solution);
    allocVecSol(&solution);
    initSolution(&solution);

    /* writting output file */
    writeSolution(&solution,&mesh);
}

void writeSolution(results * solution, grid * mesh){
    /* Auxiliary function to write the
    *  output file with solution */
    FILE *fp;

    fp = fopen(filename,"w+");
    for (int i = 0; i < imax; i++){
        fprintf(fp,"%f\t%f\t%f\t\n",mesh->x[i],solution->press[i],solution->rho[i]);
    }
    fclose(fp);
};

