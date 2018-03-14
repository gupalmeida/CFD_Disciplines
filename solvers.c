#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "input.h"
#include "mesh.h"
#include "aux.h"
#include "solvers.h"

void explicitBeamWarming(results * solution){
    /* BEAM and WARMING method using explicit
    *  Euler as time marching method and a
    *  simple central difference scheme for
    *  spatial discretization.
    *  =======================================
    *  Qnp1 = Qnp - 0.5*(dt/dx)*(Ejp1 - Ejm1)
    *  ======================================= */
    
    double Q1np1[imax], Q2np1[imax], Q3np1[imax];
    int it = 0;
    
    while (it <= itmax){
        it++;
        for (int j = 0; j < imax; j++){
            Q1np1[j] = solution->Q[0][j] - 0.5*(dt/dx)*(solution->E[0][j+1] - solution->E[0][j-1]);
            Q2np1[j] = solution->Q[1][j] - 0.5*(dt/dx)*(solution->E[1][j+1] - solution->E[1][j-1]);
            Q3np1[j] = solution->Q[2][j] - 0.5*(dt/dx)*(solution->E[2][j+1] - solution->E[2][j-1]);
        }
    }
    
}

