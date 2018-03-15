#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "input.h"
#include "aux.h"
#include "mesh.h"

void allocSolution(results * solution){
    /* Here you have memory allocation for the
    *  solution variables. The elements in vectors
    *  Q and E are accessed using the slicing
    *  in the format Q[i][j], where i refers to
    *  the equation and j to the point in the 
    *  computational mesh. For the jacobian matrix
    *  A[i][j] i and j refer to the element i, j
    *  of the jacobian matrix only. */
    solution->press = malloc(imax * sizeof(double));
    solution->vel = malloc(imax * sizeof(double));
    solution->rho = malloc(imax * sizeof(double));
    solution->energy = malloc(imax * sizeof(double));
    solution->iEnergy = malloc(imax * sizeof(double));
    solution->Q1 = malloc(imax * sizeof(double));
    solution->Q2 = malloc(imax * sizeof(double));
    solution->Q3 = malloc(imax * sizeof(double));
    solution->E1 = malloc(imax * sizeof(double));
    solution->E2 = malloc(imax * sizeof(double));
    solution->E3 = malloc(imax * sizeof(double));

    /*
    solution->A = malloc(3 * sizeof(double *));
    for (int j = 0; j < 3; j++){
        solution->A[j] = malloc(3 * sizeof(double));
    }
    */
}

void freeSolution(results * solution){
    free(solution->press);
    free(solution->vel);
    free(solution->rho);
    free(solution->energy);
    free(solution->iEnergy);
    free(solution->Q1);
    free(solution->Q2);
    free(solution->Q3);
    free(solution->E1);
    free(solution->E2);
    free(solution->E3);

    /*
    for (int j = 0; j < 3; j++){
        free(solution->A[j]);
    }
    free(solution->A);
    */
}

void initSolution(results * solution){
    int mid = (imax - 1)/2;

    for (int i = 0; i < imax; i++){
        if (i <= mid){
            solution->press[i] = p4;
            solution->rho[i] = p4;
            solution->energy[i] = 1.0/(gamma - 1.0);
            solution->iEnergy[i] = (1.0/(gamma - 1.0))/solution->rho[i];
        }
        else {
        }
    }

    for (int i = 0; i < imax; i++){
        solution->vel[i] = U0;
    }

    for (int i = 0; i < imax; i++){
        solution->Q1[i] = solution->rho[i];
        solution->Q2[i] = solution->rho[i] * solution->vel[i];
        solution->Q3[i] = solution->energy[i];
    }

    for (int i = 0; i < imax; i++){
        solution->E1[i] = solution->Q2[i];
        solution->E2[i] = (solution->Q2[i] * solution->Q2[i])/solution->Q1[i] + solution->press[i];
        solution->E3[i] = (solution->Q3[i] + solution->press[i])*solution->vel[i];
    }
}

/*
void jacobian(results * solution){
    * The jacobian matrix is determined
    *  by deriving the flux variables in
    *  relation to the conserved ones. *
    solution->A[0][0] = 0.0;
    solution->A[0][1] = 1.0;
    solution->A[0][2] = 0.0;
}
*/

void calcPrimitives(results * solution){
    double temp=0.0;
    for (int j = 0; j < imax; j++){
        solution->rho[j] = solution->Q1[j];
        solution->vel[j] = solution->Q2[j] / solution->Q1[j];
        temp = solution->Q3[j] - 0.5*(solution->Q2[j] * solution->Q2[j])/solution->Q1[j];
        solution->press[j] = temp*(gamma-1.0);
        solution->iEnergy[j] = temp/solution->Q1[j];
    }
}

void calcFluxes(results * solution){
    for (int j = 0; j < imax; j++){
        solution->E1[j] = solution->Q2[j];
        solution->E2[j] = (solution->Q2[j] * solution->Q2[j])/solution->Q1[j] + solution->press[j];
        solution->E3[j] = (solution->Q3[j] + solution->press[j])*solution->vel[j];
    }
}
