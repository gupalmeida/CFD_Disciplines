#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "input.h"
#include "mesh.h"
#include "aux.h"
#include "solvers.h"

void explicitBeamWarming(results * solution, grid * mesh){
    /* BEAM and WARMING method using explicit
    *  Euler as time marching method and a
    *  simple central difference scheme for
    *  spatial discretization.
    *  =======================================
    *  Qnp1 = Qnp - 0.5*(dt/dx)*(Ejp1 - Ejm1)
    *  ======================================= */
    
    double Q1np1[imax], Q2np1[imax], Q3np1[imax];
    int it = 0;
    double dx = mesh->x[1] - mesh->x[0];
    double dt = dx*cfl;
    
    while (it <= itmax){
        it++;
        for (int j = 0; j < imax; j++){
            Q1np1[j] = solution->Q1[j] - 0.5*(dt/dx)*(solution->E1[j+1] - solution->E1[j-1]);
            Q2np1[j] = solution->Q2[j] - 0.5*(dt/dx)*(solution->E2[j+1] - solution->E2[j-1]);
            Q3np1[j] = solution->Q3[j] - 0.5*(dt/dx)*(solution->E3[j+1] - solution->E3[j-1]);
        }
        for (int j = 0; j < imax; j++){
            solution->Q1[j] = Q1np1[j];
            solution->Q2[j] = Q2np1[j];
            solution->Q3[j] = Q3np1[j];
        }

        calcPrimitives(solution);
        calcFluxes(solution);
    }
    
}
