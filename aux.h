#ifndef AUX_H
#define AUX_H

typedef struct results{
    double *press;
    double *vel;
    double *rho;
    double *energy;
    double *iEnergy;
    double **Q;
    double **E;
    double **A;
}results;

typedef struct exact{
    double pr;
    double p2;
    double p3;
    double pl;
    double p5;
    double u1;
    double u2;
    double u3;
    double u4;
    double u5;
    double r1;
    double r2;
    double r3;
    double r4;
    double r5;
    double e1;
    double e2;
    double e3;
    double e4;
    double e5;
}exact;
/* AUXILIARY FUNCTIONS AND DEFINITIONS */

void allocSolution(results *);
void allocVecSol(results *);
void initSolution(results *);
void jacobian(results *);
void calcPrimitives(results *);
void calcFluxes(results *);

#endif
