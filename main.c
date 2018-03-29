#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "input.h"
#include "mesh.h"
#include "aux.h"
#include "solvers.h"

void writeSolution(results * , grid * ,char fileName[]);

int main(){
    
    /* parameter definition */
    grid mesh;
    results solution;
    results exact;
    
    /* generating mesh */
    meshAlloc(&mesh);
    createMesh(&mesh);
    
    /* initializing variables */
    allocSolution(&solution);
    initSolution(&solution);
    //writeSolution(&solution,&mesh);

    /* computing the exact solution */
    exactSolution(&exact,&mesh);
    writeSolution(&exact,&mesh,"exact.dat");

    /* solving Euler equations */
    switch (method){
        case 0:
            centeredScheme(&solution,&mesh);
            break;
        case 1:
            laxWendroff(&solution,&mesh);
            break;
        case 2:
            macCormack(&solution,&mesh);
            break;
        case 3:
            stegerWarming(&solution,&mesh);
            break;
        default:
            centeredScheme(&solution,&mesh);
            break;
    } 

    /* writting output file */
    writeSolution(&solution,&mesh,filename);
    freeMesh(&mesh);
    freeSolution(&solution);
    //freeSolution(&exact);
}

void writeSolution(results * solution, grid * mesh, char fileName[]){
    /* Auxiliary function to write the
    *  output file with solution */
    FILE *fp;

    fp = fopen(fileName,"w");
    if (fp != NULL){
        for (int i = 0; i < imax; i++){
            fprintf(fp,"%f\t%f\t%f\t%f\n",mesh->x[i],solution->press[i],solution->rho[i],solution->vel[i]);
        }
        fclose(fp);
    }
    else{
        printf("Could not handle the specified file\n\n");
    }
};

