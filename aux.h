#ifndef AUX_H
#define AUX_H

typedef struct results{
    double *press;
    double *vel;
    double *rho;
    double *energy;
    double *iEnergy;
    double *Q1;
    double *Q2;
    double *Q3;
    double *E1;
    double *E2;
    double *E3;
    /*
    double **A;
    */
}results;

/* AUXILIARY FUNCTIONS AND DEFINITIONS */

void allocSolution(results *);
void freeSolution(results *);
void initSolution(results *);
/*
void jacobian(results *);
*/
void calcPrimitives(results *);
void calcFluxes(results *);

#endif
