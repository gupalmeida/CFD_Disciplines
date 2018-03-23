#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "input.h"
#include "aux.h"
#include "mesh.h"

void allocSolution(results * solution){
    solution->press = malloc((double) imax * sizeof(double));
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
    solution->dissip1 = malloc(imax * sizeof(double));
    solution->dissip2 = malloc(imax * sizeof(double));
    solution->dissip3 = malloc(imax * sizeof(double));
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
    free(solution->dissip1);
    free(solution->dissip2);
    free(solution->dissip3);
}

void initSolution(results * solution){
    int mid = (imax - 1)/2;

    for (int i = 0; i < imax; i++){
        if (i <= mid){
            solution->press[i] = p4;
            solution->rho[i] = p4;
            solution->energy[i] = p4/(gamma - 1.0);
            solution->iEnergy[i] = (p4/(gamma - 1.0))/solution->rho[i];
        }
        else {
            solution->press[i] = p1;
            solution->rho[i] = p1;
            solution->energy[i] = p1/(gamma - 1.0);
            solution->iEnergy[i] = (p1/(gamma - 1.0))/solution->rho[i];
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

    return;
}

void calcPrimitives(results * solution){
    double temp=0.0;
    for (int j = 0; j < imax; j++){
        solution->rho[j] = solution->Q1[j];
        solution->vel[j] = solution->Q2[j] / solution->Q1[j];
        temp = solution->Q3[j] - 0.5*(solution->Q2[j] * solution->Q2[j])/solution->Q1[j];
        solution->press[j] = temp*(gamma-1.0);
        solution->iEnergy[j] = temp/solution->Q1[j];
    }

    return;
}

void calcFluxes(results * solution){
    for (int j = 0; j < imax; j++){
        solution->E1[j] = solution->Q2[j];
        solution->E2[j] = (solution->Q2[j] * solution->Q2[j])/solution->Q1[j] + solution->press[j];
        solution->E3[j] = (solution->Q3[j] + solution->press[j])*solution->vel[j];
    }

    return;
}

void calcDissipation(results * solution){
    switch (dissipModel){
        case 0:
            for (int j = 0; j < imax; j++){
                solution->dissip1[j] = (alpha/8.0)*(solution->Q1[j+1] - 2.0*solution->Q1[j] + solution->Q1[j-1]);
                solution->dissip2[j] = (alpha/8.0)*(solution->Q2[j+1] - 2.0*solution->Q2[j] + solution->Q2[j-1]);
                solution->dissip3[j] = (alpha/8.0)*(solution->Q3[j+1] - 2.0*solution->Q3[j] + solution->Q3[j-1]);
            }
            break;
        case 1:
            for (int j = 0; j < imax; j++){
                solution->dissip1[j] = -(alpha/8.0)*(solution->Q1[j+2] - 4.0*solution->Q1[j+1] + 6.0*solution->Q1[j] - 4.0*solution->Q1[j-1] + solution->Q1[j-2]);
                solution->dissip2[j] = -(alpha/8.0)*(solution->Q2[j+2] - 4.0*solution->Q2[j+1] + 6.0*solution->Q2[j] - 4.0*solution->Q2[j-1] + solution->Q2[j-2]);
                solution->dissip3[j] = -(alpha/8.0)*(solution->Q3[j+2] - 4.0*solution->Q3[j+1] + 6.0*solution->Q3[j] - 4.0*solution->Q3[j-1] + solution->Q3[j-2]);
            }
            break;
        default:
            for (int j = 0; j < imax; j++){
                solution->dissip1[j] = (alpha/8.0)*(solution->Q1[j+1] - 2.0*solution->Q1[j] + solution->Q1[j-1]);
                solution->dissip2[j] = (alpha/8.0)*(solution->Q2[j+1] - 2.0*solution->Q2[j] + solution->Q2[j-1]);
                solution->dissip3[j] = (alpha/8.0)*(solution->Q3[j+1] - 2.0*solution->Q3[j] + solution->Q3[j-1]);
            }
    }
}
