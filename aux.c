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
    solution->E1p = malloc(imax * sizeof(double));
    solution->E2p = malloc(imax * sizeof(double));
    solution->E3p = malloc(imax * sizeof(double));
    solution->E1m = malloc(imax * sizeof(double));
    solution->E2m = malloc(imax * sizeof(double));
    solution->E3m = malloc(imax * sizeof(double));
    solution->dissip1 = malloc(imax * sizeof(double));
    solution->dissip2 = malloc(imax * sizeof(double));
    solution->dissip3 = malloc(imax * sizeof(double));
    solution->A = malloc(3 * sizeof(double *));

    for (int j = 0; j < 3; j++){
        solution->A[j] = malloc(3 * sizeof(double));
    }
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
    free(solution->E1p);
    free(solution->E2p);
    free(solution->E3p);
    free(solution->E1m);
    free(solution->E2m);
    free(solution->E3m);
    free(solution->dissip1);
    free(solution->dissip2);
    free(solution->dissip3);
    
    for (int j = 0; j < 3; j++){
        free(solution->A[j]);
    }
    free(solution->A);
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

    for (int i = 0; i < imax; i++){
        solution->dissip1[i] = 0.0;
        solution->dissip2[i] = 0.0;
        solution->dissip3[i] = 0.0;
    }

    calcStegerFluxes(solution);

    /* initializing the jacobian matrix */
    solution->A[0][0] = 0.0;
    solution->A[0][1] = 1.0;
    solution->A[0][2] = 0.0;
    solution->A[1][0] = 0.0;
    solution->A[1][1] = 0.0;
    solution->A[1][2] = 0.0;
    solution->A[2][0] = 0.0;
    solution->A[2][1] = 0.0;
    solution->A[2][2] = 0.0;
    
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

void calcStegerFluxes(results * solution){
    double a, u, p, rho;
    double lbd1, lbd2, lbd3;
    double lbd1p, lbd2p, lbd3p;
    double lbd1m, lbd2m, lbd3m;
    double gm = gamma;

    for (int j = 0; j < imax; j++){
        rho = solution->rho[j];
        u = solution->vel[j];
        p = solution->press[j];

        a = sqrt(gamma*(p/rho));

        lbd1 = u;
        lbd2 = u + a;
        lbd3 = u - a;

        /* computing E+ fluxes */
        lbd1p = 0.5*(lbd1 + fabs(lbd1));
        lbd2p = 0.5*(lbd2 + fabs(lbd2));
        lbd3p = 0.5*(lbd3 + fabs(lbd3));

        solution->E1p[j] = (rho/(2.0*gm))*(2.0*(gm-1.0)*lbd1p
                + lbd2p + lbd3p);
        solution->E2p[j] = (rho/(2.0*gm))*(2.0*(gm-1.0)*lbd1p*u
                + lbd2p*(u+a) + lbd3p*(u-a));
        solution->E3p[j] = (rho/(2.0*gm))*((gm-1.0)*lbd1p*u*u
                + 0.5*lbd2p*(u+a)*(u+a)
                + 0.5*lbd3p*(u-a)*(u-a)
                + ((3.0-gm)*(lbd2p + lbd3p)*a*a)/(2.0*(gm-1.0)));

        /* computing E- fluxes */
        lbd1m = 0.5*(lbd1 - fabs(lbd1));
        lbd2m = 0.5*(lbd2 - fabs(lbd2));
        lbd3m = 0.5*(lbd3 - fabs(lbd3));

        solution->E1m[j] = (rho/(2.0*gm))*(2.0*(gm-1.0)*lbd1m
                + lbd2m + lbd3m);
        solution->E2m[j] = (rho/(2.0*gm))*(2.0*(gm-1.0)*lbd1m*u
                + lbd2m*(u+a) + lbd3m*(u-a));
        solution->E3m[j] = (rho/(2.0*gm))*((gm-1.0)*lbd1m*u*u
                + 0.5*lbd2m*(u+a)*(u+a)
                + 0.5*lbd3m*(u-a)*(u-a)
                + ((3.0-gm)*(lbd2m + lbd3m)*a*a)/(2.0*(gm-1.0)));
    }
    //printf("l1: %f, l2: %f, l3: %f\n",lbd1,lbd2,lbd3);

    return;
}

void calcJacobian(double q1, double q2, double q3, results * solution){
    double gm = gamma;

    solution->A[0][0] = 0.0;
    solution->A[0][1] = 1.0;
    solution->A[0][2] = 0.0;
    solution->A[1][0] = (gm - 3.0)*(q2 * q2)/(2.0 * q1 * q1);
    solution->A[1][1] = (3.0 - gm)*(q2 / q1);
    solution->A[1][2] = (gm - 1.0);
    solution->A[2][0] = (gm - 1.0)*(q2*q2*q2)/(q1*q1*q1) - gm*(q3*q2)/(q1*q1);
    solution->A[2][1] = gm*(q3/q1) - 1.5*(gm - 1.0)*(q2*q2)/(q1*q1);
    solution->A[2][2] = gm*(q2/q1);
} 

void calcDissipation(results * solution, double lambda){
    /* case 0 - 2nd order linear artificial dissipation
       case 1 - 4th order linear artificial dissipation
       case 2 - Jameson non-linear artificial dissipation
       default - 2nd order linear artificial dissipation */

    double nu[imax];
    double k2, k4;
    double dph, dmh;
    k2 = 1.0/2.0;
    k4 = 1.0/20.0;

    for (int j = 0; j < imax; j++){
        nu[j] = 0.0;
    }

    switch (dissipModel){
        case 0:
            for (int j = 2; j < imax-2; j++){
                solution->dissip1[j] = (alpha/8.0)*(1.0/lambda)*(solution->Q1[j+1] - 2.0*solution->Q1[j] + solution->Q1[j-1]);
                solution->dissip2[j] = (alpha/8.0)*(1.0/lambda)*(solution->Q2[j+1] - 2.0*solution->Q2[j] + solution->Q2[j-1]);
                solution->dissip3[j] = (alpha/8.0)*(1.0/lambda)*(solution->Q3[j+1] - 2.0*solution->Q3[j] + solution->Q3[j-1]);
            }
            break;

        case 1:

            for (int j = 2; j < imax-2; j++){
                solution->dissip1[j] = -(alpha/8.0)*(1.0/lambda)*(solution->Q1[j+2] - 4.0*solution->Q1[j+1] + 6.0*solution->Q1[j] - 4.0*solution->Q1[j-1] + solution->Q1[j-2]);
                solution->dissip2[j] = -(alpha/8.0)*(1.0/lambda)*(solution->Q2[j+2] - 4.0*solution->Q2[j+1] + 6.0*solution->Q2[j] - 4.0*solution->Q2[j-1] + solution->Q2[j-2]);
                solution->dissip3[j] = -(alpha/8.0)*(1.0/lambda)*(solution->Q3[j+2] - 4.0*solution->Q3[j+1] + 6.0*solution->Q3[j] - 4.0*solution->Q3[j-1] + solution->Q3[j-2]);
            }

            break;

        case 2:

            for (int j = 2; j < imax-2; j++){
                nu[j] = fabs(solution->press[j+1] - 2.0*solution->press[j] + solution->press[j-1]) / (fabs(solution->press[j+1]) + 2.0*fabs(solution->press[j]) + fabs(solution->press[j-1]));
            }

            for (int j = 2; j < imax-2; j++){
                double eps2ph = k2*max(nu[j+1],nu[j]);
                double eps4ph = max(0.0, (k4 - eps2ph) );
                double eps2mh = k2*max(nu[j],nu[j-1]);
                double eps4mh = max(0.0, (k4 - eps2mh) );
                dph = eps2ph*(solution->Q1[j+1] - solution->Q1[j]) - eps4ph*(solution->Q1[j+2] - 3.0*solution->Q1[j+1] + 3.0*solution->Q1[j] - solution->Q1[j-1]);
                dmh = eps2mh*(solution->Q1[j] - solution->Q1[j-1]) - eps4mh*(solution->Q1[j+1] - 3.0*solution->Q1[j] + 3.0*solution->Q1[j-1] - solution->Q1[j-2]);
                solution->dissip1[j] = (1.0/lambda)*(dph - dmh);
            }

            for (int j = 2; j < imax-2; j++){
                double eps2ph = k2*max(nu[j+1],nu[j]);
                double eps4ph = max(0.0, (k4 - eps2ph) );
                double eps2mh = k2*max(nu[j],nu[j-1]);
                double eps4mh = max(0.0, (k4 - eps2mh) );
                dph = eps2ph*(solution->Q2[j+1] - solution->Q2[j]) - eps4ph*(solution->Q2[j+2] - 3.0*solution->Q2[j+1] + 3.0*solution->Q2[j] - solution->Q2[j-1]);
                dmh = eps2mh*(solution->Q2[j] - solution->Q2[j-1]) - eps4mh*(solution->Q2[j+1] - 3.0*solution->Q2[j] + 3.0*solution->Q2[j-1] - solution->Q2[j-2]);
                solution->dissip2[j] = (1.0/lambda)*(dph - dmh);
            }
            
            for (int j = 2; j < imax-2; j++){
                double eps2ph = k2*max(nu[j+1],nu[j]);
                double eps4ph = max(0.0, (k4 - eps2ph) );
                double eps2mh = k2*max(nu[j],nu[j-1]);
                double eps4mh = max(0.0, (k4 - eps2mh) );
                dph = eps2ph*(solution->Q3[j+1] - solution->Q3[j]) - eps4ph*(solution->Q3[j+2] - 3.0*solution->Q3[j+1] + 3.0*solution->Q3[j] - solution->Q3[j-1]);
                dmh = eps2mh*(solution->Q3[j] - solution->Q3[j-1]) - eps4mh*(solution->Q3[j+1] - 3.0*solution->Q3[j] + 3.0*solution->Q3[j-1] - solution->Q3[j-2]);
                solution->dissip3[j] = (1.0/lambda)*(dph - dmh);
            }

            break;

        default:
            for (int j = 2; j < imax-2; j++){
                solution->dissip1[j] = (alpha/8.0)*(1.0/lambda)*(solution->Q1[j+1] - 2.0*solution->Q1[j] + solution->Q1[j-1]);
                solution->dissip2[j] = (alpha/8.0)*(1.0/lambda)*(solution->Q2[j+1] - 2.0*solution->Q2[j] + solution->Q2[j-1]);
                solution->dissip3[j] = (alpha/8.0)*(1.0/lambda)*(solution->Q3[j+1] - 2.0*solution->Q3[j] + solution->Q3[j-1]);
            }
            break;
    }

	//free(nu);
	//free(eps2);
	//free(eps4);
}

double max(double a, double b){
    if (a > b){
        return a;
    }
    else{
        return b;
    }
}
