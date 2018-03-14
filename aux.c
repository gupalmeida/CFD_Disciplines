#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "input.h"
#include "aux.h"
#include "mesh.h"

void allocSolution(results * solution){
    solution->press = malloc(imax * sizeof(double));
    solution->vel = malloc(imax * sizeof(double));
    solution->rho = malloc(imax * sizeof(double));
    solution->energy = malloc(imax * sizeof(double));
    solution->iEnergy = malloc(imax * sizeof(double));
}

void allocVecSol(results * solution){
    /* Here you have memory allocation for the
    *  solution vectors Q (conserved variables)
    *  and E (flux variables). The elements are
    *  accessed using the slicing in the format
    *  Q[i][j], where i refers to the equation
    *  and j to the point in the computational
    *  mesh. For the jacobian matrix A[i][j]
    *  i and j refer to the element i, j of the
    *  jacobian matrix only. */
    solution->Q = malloc(3 * sizeof(double *));
    for (int j = 0; j < 3; j++){
        solution->Q[j] = malloc(imax * sizeof(double));
    }

    solution->E = malloc(3 * sizeof(double *));
    for (int j = 0; j < 3; j++){
        solution->E[j] = malloc(imax * sizeof(double));
    }

    solution->A = malloc(3 * sizeof(double *));
    for (int j = 0; j < 3; j++){
        solution->E[j] = malloc(3 * sizeof(double));
    }
}

void initSolution(results * solution){
    double mid = (double) (imax - 1)/2;

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
}

void jacobian(results * solution){
    /* The jacobian matrix is determined
    *  by deriving the flux variables in
    *  relation to the conserved ones. */
    solution->A[0][0] = 0.0;
    solution->A[0][1] = 1.0;
    solution->A[0][2] = 0.0;
}

void calcPrimitives(results * solution){
    double temp=0.0;
    for (int j = 0; j < imax; j++){
        solution->rho[j] = solution->Q[0][j];
        solution->vel[j] = pow(solution->Q[1][j],2.0) / solution->Q[0][j];
        temp = solution->Q[1][j] - 0.5*pow(solution->Q[1][j],2.0)/solution->Q[0][j];
        solution->press[j] = temp*(gamma-1.0);
        solution->iEnergy[j] = temp/solution->Q[0][j];
    }
}

void calcFluxes(results * solution){
    for (int j = 0; j < imax; j++){
        solution->E[0][j] = solution->Q[1][j];
        solution->E[1][j] = pow(solution->Q[1][j],2.0)/solution->Q[0][j] + solution->press[j];
        solution->E[2][j] = (solution->Q[2][j] + solution->press[j])*solution->vel[j];
    }
}
