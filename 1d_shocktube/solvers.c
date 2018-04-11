#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "input.h"
#include "mesh.h"
#include "aux.h"
#include "solvers.h"

/* beam and warming explicit centered scheme */
void centeredScheme(results * solution, grid * mesh){
    /* BEAM and WARMING method using explicit
    *  Euler as time marching method and a
    *  simple central difference scheme for
    *  spatial discretization.
    *  =======================================
    *  Qnp1 = Qnp - 0.5*(dt/dx)*(Ejp1 - Ejm1)
    *  ======================================= */
    
    double *Q1np1 = malloc(imax * sizeof(double));
    double *Q2np1 = malloc(imax * sizeof(double));
    double *Q3np1 = malloc(imax * sizeof(double));
    int it = 0;
    double dx = mesh->x[1] - mesh->x[0];
    double a = sqrt(gamma);
    double dt = (dx*cfl)/a;
    double t = 0.0;
    double lambda = dt/dx;

    printf("\n\n===============================\n");
    printf("         Centered Scheme\n");
    printf("===============================\n\n");
    
    while (t <= (double) tmax && it < (int) itmax){
        /* calculates the flux vector */
        calcFluxes(solution);

        /* calculates dissipation for the current solution Q */
        calcDissipation(solution, lambda);

        /* marching the solution */
        for (int j = 2; j < imax-2; j++){
            Q1np1[j] = solution->Q1[j] - 0.5*lambda*(solution->E1[j+1] - solution->E1[j-1]) + lambda*(solution->dissip1[j]);
            Q2np1[j] = solution->Q2[j] - 0.5*lambda*(solution->E2[j+1] - solution->E2[j-1]) + lambda*(solution->dissip2[j]);
            Q3np1[j] = solution->Q3[j] - 0.5*lambda*(solution->E3[j+1] - solution->E3[j-1]) + lambda*(solution->dissip3[j]);
        }
        for (int j = 2; j < imax-2; j++){
            solution->Q1[j] = Q1np1[j];
            solution->Q2[j] = Q2np1[j];
            solution->Q3[j] = Q3np1[j];
        }

        calcPrimitives(solution);
        t = t + dt;
        it++;
    
        if ( (it%printAt == 0 ) || (fmod(tmax,t) >= 1.0)){
            printf("Iteration: %d, simulation time: %lf\n",it,t);
        }
    }

    /* deallocating variables */
    free(Q1np1);
    free(Q2np1);
    free(Q3np1);

    printf("\n\n");

}

/* lax-wendroff method */
void laxWendroff(results * solution, grid * mesh){
    /* LAX-WENDROFF method using explicit
    *  Euler as time marching method.
    *  =======================================
    *  Qjnp1 = Qjn - 0.5*lambda*(Ejp1 - Ejm1)
    *        + 0.5*lambda2*(Ajph*(Ejp1 - Ej))
    *        - 0.5*lambda2*(Ajmh*(Ej - Ejm1))
    *  ======================================= */
    
    double *Q1np1 = malloc(imax * sizeof(double));
    double *Q2np1 = malloc(imax * sizeof(double));
    double *Q3np1 = malloc(imax * sizeof(double));
    int it = 0;
    double dx = mesh->x[1] - mesh->x[0];
    double q1, q2, q3;
    double temp1, temp2, temp3;
    double a = sqrt(gamma);
    double dt = (dx*cfl)/a;
    double t = 0.0;
    double lambda = dt/dx;

    printf("\n\n===============================\n");
    printf("       Lax-Wendroff Scheme\n");
    printf("===============================\n\n");
    
    while (t <= (double) tmax && it < (int) itmax){
        /* calculates the flux vector */
        calcFluxes(solution);

        /* calculates dissipation for the current solution Q */
        calcDissipation(solution, lambda);

        /* marching the solution */
        for (int j = 2; j < imax-2; j++){
            /* calculating the variables in point j+1/2 */
            q1 = 0.5*(solution->Q1[j] + solution->Q1[j+1]);
            q2 = 0.5*(solution->Q2[j] + solution->Q2[j+1]);
            q3 = 0.5*(solution->Q3[j] + solution->Q3[j+1]);

            calcJacobian(q1,q2,q3,solution);

            temp1 = solution->Q1[j] 
                - 0.5*lambda*(solution->E1[j+1] - solution->E1[j-1]) 
                + 0.5*lambda*lambda*(solution->A[0][0]*(solution->E1[j+1]-solution->E1[j]) + solution->A[0][1]*(solution->E2[j+1]-solution->E2[j]) + solution->A[0][2]*(solution->E3[j+1]-solution->E3[j]));

            temp2 = solution->Q2[j] 
                - 0.5*lambda*(solution->E2[j+1] - solution->E2[j-1]) 
                + 0.5*lambda*lambda*(solution->A[1][0]*(solution->E1[j+1]-solution->E1[j]) + solution->A[1][1]*(solution->E2[j+1]-solution->E2[j]) + solution->A[1][2]*(solution->E3[j+1]-solution->E3[j]));

            temp3 = solution->Q3[j] 
                - 0.5*lambda*(solution->E3[j+1] - solution->E3[j-1]) 
                + 0.5*lambda*lambda*(solution->A[2][0]*(solution->E1[j+1]-solution->E1[j]) + solution->A[2][1]*(solution->E2[j+1]-solution->E2[j]) + solution->A[2][2]*(solution->E3[j+1]-solution->E3[j]));

            /* calculating the variables in point j-1/2 */
            q1 = 0.5*(solution->Q1[j] + solution->Q1[j-1]);
            q2 = 0.5*(solution->Q2[j] + solution->Q2[j-1]);
            q3 = 0.5*(solution->Q3[j] + solution->Q3[j-1]);

            calcJacobian(q1,q2,q3,solution);

            Q1np1[j] = temp1
                - 0.5*lambda*lambda*(solution->A[0][0]*(solution->E1[j]-solution->E1[j-1]) + solution->A[0][1]*(solution->E2[j]-solution->E2[j-1]) + solution->A[0][2]*(solution->E3[j]-solution->E3[j-1]))
                + lambda*(solution->dissip1[j]);

            Q2np1[j] = temp2
                - 0.5*lambda*lambda*(solution->A[1][0]*(solution->E1[j]-solution->E1[j-1]) + solution->A[1][1]*(solution->E2[j]-solution->E2[j-1]) + solution->A[1][2]*(solution->E3[j]-solution->E3[j-1]))
                + lambda*(solution->dissip2[j]);

            Q3np1[j] = temp3
                - 0.5*lambda*lambda*(solution->A[2][0]*(solution->E1[j]-solution->E1[j-1]) + solution->A[2][1]*(solution->E2[j]-solution->E2[j-1]) + solution->A[2][2]*(solution->E3[j]-solution->E3[j-1]))
                + lambda*(solution->dissip3[j]);
        }
        for (int j = 2; j < imax-2; j++){
            solution->Q1[j] = Q1np1[j];
            solution->Q2[j] = Q2np1[j];
            solution->Q3[j] = Q3np1[j];
        }

        calcPrimitives(solution);
        t = t + dt;
        it++;
    
        if ( (it%printAt == 0 ) || (fmod(tmax,t) >= 1.0)){
            printf("Iteration: %d, simulation time: %lf\n",it,t);
        }
    }

    /* deallocating variables */
    free(Q1np1);
    free(Q2np1);
    free(Q3np1);

    printf("\n\n");

}

/* classical predictor-corrector macCormack scheme */
void macCormack(results * solution, grid * mesh){
    /* MACCORMACK method using explicit
    *  Euler as time marching method.
    *  ==============================================
    *  predictor step
    *  bQjnp1 = Qjn - lambda*(Ejp1 - Ej)
    *  
    *  corrector step
    *  Qjnp1 = 0.5*(Qjn + bQjnp1 - lambda*(Ej - Ejm1)
    *  ============================================== */
    
    double *Q1np1 = malloc(imax * sizeof(double));
    double *Q2np1 = malloc(imax * sizeof(double));
    double *Q3np1 = malloc(imax * sizeof(double));
    double bQ1np1, bQ2np1, bQ3np1;
    int it = 0;
    double dx = mesh->x[1] - mesh->x[0];
    double a = sqrt(gamma);
    double dt = (dx*cfl)/a;
    double t = 0.0;
    double lambda = dt/dx;

    printf("\n\n===============================\n");
    printf("         MacCormack Scheme\n");
    printf("===============================\n\n");
    
    while (t <= (double) tmax && it < (int) itmax){
        /* calculates the flux vector */
        calcFluxes(solution);

        /* calculates dissipation for the current solution Q */
        calcDissipation(solution, lambda);

        /* marching the solution */
        for (int j = 2; j < imax-2; j++){
            /* predictor step */
            bQ1np1 = solution->Q1[j] - lambda*(solution->E1[j+1] - solution->E1[j]);
            bQ2np1 = solution->Q2[j] - lambda*(solution->E2[j+1] - solution->E2[j]);
            bQ3np1 = solution->Q3[j] - lambda*(solution->E3[j+1] - solution->E3[j]);

            /* corrector step */
            Q1np1[j] = 0.5*(solution->Q1[j] + bQ1np1 - lambda*(solution->E1[j] - solution->E1[j-1])) + lambda*(solution->dissip1[j]);
            Q2np1[j] = 0.5*(solution->Q2[j] + bQ2np1 - lambda*(solution->E2[j] - solution->E2[j-1])) + lambda*(solution->dissip2[j]);
            Q3np1[j] = 0.5*(solution->Q3[j] + bQ3np1 - lambda*(solution->E3[j] - solution->E3[j-1])) + lambda*(solution->dissip3[j]);
        }
        for (int j = 2; j < imax-2; j++){
            solution->Q1[j] = Q1np1[j];
            solution->Q2[j] = Q2np1[j];
            solution->Q3[j] = Q3np1[j];
        }

        calcPrimitives(solution);
        t = t + dt;
        it++;
    
        if ( (it%printAt == 0 ) || (fmod(tmax,t) >= 1.0)){
            printf("Iteration: %d, simulation time: %lf\n",it,t);
        }
    }

    /* deallocating variables */
    free(Q1np1);
    free(Q2np1);
    free(Q3np1);

    printf("\n\n");

}

/* steger and warming flux vector splitting scheme */
void stegerWarming(results * solution, grid * mesh){
    /* STEGER AND WARMING flux vector splitting scheme
    *  using explicit Euler as time marching method.
    *  ==============================================
    *  Qjnp1 = Qjn -
    *  ============================================== */
    
    double *Q1np1 = malloc(imax * sizeof(double));
    double *Q2np1 = malloc(imax * sizeof(double));
    double *Q3np1 = malloc(imax * sizeof(double));
    int it = 0;
    double dx = mesh->x[1] - mesh->x[0];
    double a = sqrt(gamma);
    double dt = (dx*cfl)/a;
    double t = 0.0;
    double lambda = dt/dx;

    printf("\n\n=================================\n");
    printf("  Steger and Warming FVS Scheme\n");
    printf("=================================\n\n");
    
    while (t <= (double) tmax && it < (int) itmax){
        /* calculates the flux vector */
        calcStegerFluxes(solution);

        /* calculates dissipation for the current solution Q */
        //calcDissipation(solution, lambda);

        /* marching the solution */
        for (int j = 2; j < imax-2; j++){
            switch (order){
                case 1:
                    Q1np1[j] = solution->Q1[j]
                        - lambda*(solution->E1p[j] - solution->E1p[j-1])
                        - lambda*(solution->E1m[j+1] - solution->E1m[j]);
                    Q2np1[j] = solution->Q2[j]
                        - lambda*(solution->E2p[j] - solution->E2p[j-1])
                        - lambda*(solution->E2m[j+1] - solution->E2m[j]);
                    Q3np1[j] = solution->Q3[j]
                        - lambda*(solution->E3p[j] - solution->E3p[j-1])
                        - lambda*(solution->E3m[j+1] - solution->E3m[j]);
                    break;
                case 2:
                    Q1np1[j] = solution->Q1[j]
                        - (lambda/2.0)*(3.0*solution->E1p[j] - 4.0*solution->E1p[j-1] + solution->E1p[j-2])
                        - (lambda/2.0)*(-3.0*solution->E1m[j] + 4.0*solution->E1m[j+1] - solution->E1m[j+2]);
                    Q2np1[j] = solution->Q2[j]
                        - (lambda/2.0)*(3.0*solution->E2p[j] - 4.0*solution->E2p[j-1] + solution->E2p[j-2])
                        - (lambda/2.0)*(-3.0*solution->E2m[j] + 4.0*solution->E2m[j+1] - solution->E2m[j+2]);
                    Q3np1[j] = solution->Q3[j]
                        - (lambda/2.0)*(3.0*solution->E3p[j] - 4.0*solution->E3p[j-1] + solution->E3p[j-2])
                        - (lambda/2.0)*(-3.0*solution->E3m[j] + 4.0*solution->E3m[j+1] - solution->E3m[j+2]);
                    break;
                default:
                    Q1np1[j] = solution->Q1[j]
                        - lambda*(solution->E1p[j] - solution->E1p[j-1])
                        - lambda*(solution->E1m[j+1] - solution->E1m[j]);
                    Q2np1[j] = solution->Q2[j]
                        - lambda*(solution->E2p[j] - solution->E2p[j-1])
                        - lambda*(solution->E2m[j+1] - solution->E2m[j]);
                    Q3np1[j] = solution->Q3[j]
                        - lambda*(solution->E3p[j] - solution->E3p[j-1])
                        - lambda*(solution->E3m[j+1] - solution->E3m[j]);
                    break;
            }
        }

        for (int j = 2; j < imax-2; j++){
            solution->Q1[j] = Q1np1[j];
            solution->Q2[j] = Q2np1[j];
            solution->Q3[j] = Q3np1[j];
        }

        calcPrimitives(solution);
        t = t + dt;
        it++;
    
        if ( (it%printAt == 0 ) || (fmod(tmax,t) >= 1.0)){
            printf("Iteration: %d, simulation time: %lf\n",it,t);
        }
    }

    /* deallocating variables */
    free(Q1np1);
    free(Q2np1);
    free(Q3np1);

    printf("\n\n");

}

/* van Leer non-MUSCL flux vector splitting scheme */
void vanLeerNonMUSCL(results * solution, grid * mesh){
    /* STEGER AND WARMING flux vector splitting scheme
    *  using explicit Euler as time marching method.
    *  ==============================================
    *  Qjnp1 = Qjn -
    *  ============================================== */
    
    double *Q1np1 = malloc(imax * sizeof(double));
    double *Q2np1 = malloc(imax * sizeof(double));
    double *Q3np1 = malloc(imax * sizeof(double));
    int it = 0;
    double dx = mesh->x[1] - mesh->x[0];
    double a = sqrt(gamma);
    double dt = (dx*cfl)/a;
    double t = 0.0;
    double lambda = dt/dx;

    printf("\n\n=================================\n");
    printf("  van Leer non-MUSCL FVS Scheme\n");
    printf("=================================\n\n");
    
    while (t <= (double) tmax && it < (int) itmax){
        /* calculates the flux vector */
        calcVanLeerFluxes(solution);

        /* calculates dissipation for the current solution Q */
        //calcDissipation(solution, lambda);

        /* marching the solution */
        for (int j = 2; j < imax-2; j++){
            switch (order){
                case 1:
                    Q1np1[j] = solution->Q1[j]
                        - lambda*(solution->E1p[j] - solution->E1p[j-1])
                        - lambda*(solution->E1m[j+1] - solution->E1m[j]);
                    Q2np1[j] = solution->Q2[j]
                        - lambda*(solution->E2p[j] - solution->E2p[j-1])
                        - lambda*(solution->E2m[j+1] - solution->E2m[j]);
                    Q3np1[j] = solution->Q3[j]
                        - lambda*(solution->E3p[j] - solution->E3p[j-1])
                        - lambda*(solution->E3m[j+1] - solution->E3m[j]);
                    break;
                case 2:
                    Q1np1[j] = solution->Q1[j]
                        - (lambda/2.0)*(3.0*solution->E1p[j] - 4.0*solution->E1p[j-1] + solution->E1p[j-2])
                        - (lambda/2.0)*(-3.0*solution->E1m[j] + 4.0*solution->E1m[j+1] - solution->E1m[j+2]);
                    Q2np1[j] = solution->Q2[j]
                        - (lambda/2.0)*(3.0*solution->E2p[j] - 4.0*solution->E2p[j-1] + solution->E2p[j-2])
                        - (lambda/2.0)*(-3.0*solution->E2m[j] + 4.0*solution->E2m[j+1] - solution->E2m[j+2]);
                    Q3np1[j] = solution->Q3[j]
                        - (lambda/2.0)*(3.0*solution->E3p[j] - 4.0*solution->E3p[j-1] + solution->E3p[j-2])
                        - (lambda/2.0)*(-3.0*solution->E3m[j] + 4.0*solution->E3m[j+1] - solution->E3m[j+2]);
                    break;
                default:
                    Q1np1[j] = solution->Q1[j]
                        - lambda*(solution->E1p[j] - solution->E1p[j-1])
                        - lambda*(solution->E1m[j+1] - solution->E1m[j]);
                    Q2np1[j] = solution->Q2[j]
                        - lambda*(solution->E2p[j] - solution->E2p[j-1])
                        - lambda*(solution->E2m[j+1] - solution->E2m[j]);
                    Q3np1[j] = solution->Q3[j]
                        - lambda*(solution->E3p[j] - solution->E3p[j-1])
                        - lambda*(solution->E3m[j+1] - solution->E3m[j]);
                    break;
            }
        }

        for (int j = 2; j < imax-2; j++){
            solution->Q1[j] = Q1np1[j];
            solution->Q2[j] = Q2np1[j];
            solution->Q3[j] = Q3np1[j];
        }

        calcPrimitives(solution);
        t = t + dt;
        it++;
    
        if ( (it%printAt == 0 ) || (fmod(tmax,t) >= 1.0)){
            printf("Iteration: %d, simulation time: %lf\n",it,t);
        }
    }

    /* deallocating variables */
    free(Q1np1);
    free(Q2np1);
    free(Q3np1);

    printf("\n\n");

}

/* liou AUSM+ flux vector splitting scheme */
void liouAUSMplus(results * solution, grid * mesh){
    /* STEGER AND WARMING flux vector splitting scheme
    *  using explicit Euler as time marching method.
    *  ==============================================
    *  Qjnp1 = Qjn -
    *  ============================================== */
    
    double *Q1np1 = malloc(imax * sizeof(double));
    double *Q2np1 = malloc(imax * sizeof(double));
    double *Q3np1 = malloc(imax * sizeof(double));
    int it = 0;
    double dx = mesh->x[1] - mesh->x[0];
    double a = sqrt(gamma);
    double dt = (dx*cfl)/a;
    double t = 0.0;
    double lambda = dt/dx;

    printf("\n\n=================================\n");
    printf("     Liou AUSM+ FVS Scheme\n");
    printf("=================================\n\n");
    
    while (t <= (double) tmax && it < (int) itmax){
        /* calculates the flux vector */
        calcFluxes(solution);
        calcLiouFluxes(solution);

        /* calculates dissipation for the current solution Q */
        //calcDissipation(solution, lambda);

        /* marching the solution */
        for (int j = 2; j < imax-2; j++){
            Q1np1[j] = solution->Q1[j]
                - lambda*(solution->E1p[j] - solution->E1p[j-1]);
            Q2np1[j] = solution->Q2[j]
                - lambda*(solution->E2p[j] - solution->E2p[j-1]);
            Q3np1[j] = solution->Q3[j]
                - lambda*(solution->E3p[j] - solution->E3p[j-1]);
        }

        for (int j = 2; j < imax-2; j++){
            solution->Q1[j] = Q1np1[j];
            solution->Q2[j] = Q2np1[j];
            solution->Q3[j] = Q3np1[j];
        }

        calcPrimitives(solution);
        t = t + dt;
        it++;
    
        if ( (it%printAt == 0 ) || (fmod(tmax,t) >= 1.0)){
            printf("Iteration: %d, simulation time: %lf\n",it,t);
        }
    }

    /* deallocating variables */
    free(Q1np1);
    free(Q2np1);
    free(Q3np1);

    printf("\n\n");

}

/* Roe approximate Riemann solver method */
void roeMethod(results * solution, grid * mesh){
    /* ROE approximate riemann solver method
    *  using explicit Euler as time marching method.
    *  ==============================================
    *  Qjnp1 = Qjn -
    *  ============================================== */
    
    double *Q1np1 = malloc(imax * sizeof(double));
    double *Q2np1 = malloc(imax * sizeof(double));
    double *Q3np1 = malloc(imax * sizeof(double));
    int it = 0;
    double dx = mesh->x[1] - mesh->x[0];
    double a = sqrt(gamma);
    double dt = (dx*cfl)/a;
    double t = 0.0;
    double lambda = dt/dx;

    printf("\n\n==========================================\n");
    printf("     Roe\'s approximate Riemann solver     \n");
    printf("==========================================\n");
    
    while (t <= (double) tmax && it < (int) itmax){
        /* calculates the flux vector */
        calcFluxes(solution);
        calcRoeFluxes(solution);

        /* calculates dissipation for the current solution Q */
        //calcDissipation(solution, lambda);

        /* marching the solution */
        for (int j = 2; j < imax-2; j++){
            Q1np1[j] = solution->Q1[j]
                - lambda*(solution->E1p[j] - solution->E1p[j-1]);
            Q2np1[j] = solution->Q2[j]
                - lambda*(solution->E2p[j] - solution->E2p[j-1]);
            Q3np1[j] = solution->Q3[j]
                - lambda*(solution->E3p[j] - solution->E3p[j-1]);
        }

        for (int j = 2; j < imax-2; j++){
            solution->Q1[j] = Q1np1[j];
            solution->Q2[j] = Q2np1[j];
            solution->Q3[j] = Q3np1[j];
        }

        calcPrimitives(solution);
        t = t + dt;
        it++;
    
        if ( (it%printAt == 0 ) || (fmod(tmax,t) >= 1.0)){
            printf("Iteration: %d, simulation time: %lf\n",it,t);
        }
    }

    /* deallocating variables */
    free(Q1np1);
    free(Q2np1);
    free(Q3np1);

    printf("\n\n");

}

/* exact solution */
void exactSolution (results * solution, grid * mesh){

	double pr1,pr2,pr3,pr4,pr6;
	double rhor1,rhor2,rhor3,rhor4,rhor6;
	double ur1,ur2,ur3,ur4,ur6;
	double A6, B6;
	double g = gamma;
	double x, x0, t;
	double a1, a3, a6;
	double u_tail, u_head, u_contact, u_shock;
	double f, fl, pi, pOld, pNew;

	solution->rho = malloc (imax * sizeof(double));
	solution->press = malloc (imax * sizeof(double));
	solution->vel = malloc (imax * sizeof(double));

	x0 = 0.0;
	t = tmax;

	pr1 = p4;
	rhor1 = p4;
	ur1 = 0.0;
	a1 = sqrt ( g * pr1/rhor1 );

	rhor6 = p1;
	pr6 = p1;
	ur6 = 0.0;
	a6 = sqrt ( g * pr6/rhor6 );

	A6 = 2.0 / ( rhor6 * (g + 1.0) );
	B6 = pr6 * (g-1.0) / (g+1.0);

	pi = (pr1+pr6)/10.0;
	pNew = pi;
	pOld = 0.0;

	f = ( pi - pr6 ) * sqrt( A6 / (pi + B6)) + 2.0 * a1 / (g-1.0) * ( pow( (pi/pr1),((g-1.0)/(2.0*g)) ) - 1.0);

	fl = sqrt( A6 / (pi + B6) ) - 0.5 * A6 * (pi - pr6) * sqrt( (pi + B6)/A6 ) / ((pi + B6)*(pi + B6)) + a1 / (g*pr1) * (pow ((pi/pr1),( (-g-1.0)/(2.0*g) )) );

    // iterative method
	while (fabs(pNew - pOld) > tol){
		pOld = pi;
		pNew = pi - f / fl;
		pi = pNew;

		f = ( pi - pr6 ) * sqrt( A6 / (pi + B6)) + 2.0 * a1 / (g-1.0) * ( pow( (pi/pr1),((g-1.0)/(2.0*g)) ) - 1.0 ) + ur6 - ur1;

		fl = sqrt( A6 / (pi + B6) ) - 0.5 * A6 * (pi - pr6) * sqrt( (pi + B6)/A6 ) / ((pi + B6)*(pi + B6)) + a1 / (g*pr1) * (pow ((pi/pr1),( (-g-1.0)/(2.0*g) )) );
	}


	pr3 = pNew;

	rhor3 = rhor1 * pow( (pr3/pr1), 1.0/g );

	ur3 = ur1 - ( 2.0 * a1 / (g - 1.0) ) * (pow( (pr3/pr1), ((g-1.0)/(2.0*g)) ) -1.0 );

	a3 = sqrt( g * pr3 / rhor3 );

	pr4 = pr3;

	rhor4 = rhor6 * ((pr6 * (g-1.0) + pr4 * (g+1.0)) / (pr4 * (g-1.0) + pr6 * (g+1.0)) );

	ur4 = ur6 + (pr4 - pr6) * sqrt( A6 / (pr4 + B6)) ;

	u_head = ur1 - a1;

	u_tail = ur3 - a3;

	u_contact = ur3;

	u_shock = ur6 + a6 * sqrt( (g+1.0)*pr4 / (2.0 * g * pr6) + (g-1.0)/(2.0*g) );


	for ( int i = 0 ; i < imax; i++ ){
		x = mesh->x[i];

		rhor2 = rhor1 * pow( ( 2.0 /(g+1.0) + (g-1.0)/(a1*(g+1.0)) * (ur1 - (x - x0)/t) ),(2.0/(g-1.0)) );

		pr2 = pr1 * pow( ( 2.0 /(g+1.0) + (g-1.0)/(a1*(g+1.0)) * (ur1 - (x - x0)/t) ),(2.0 * g/(g-1.0)) );

		ur2 = (2.0 / (g + 1.0) ) * (a1 + 0.5 * (g-1.0) * ur1 + (x-x0)/t);

		// Region 1
		if ( ( x - x0 ) <= t * u_head){
			solution->press[i] = pr1;
			solution->rho[i] = rhor1;
			solution->vel[i] = ur1;
			//solution->s[i] = solution->press[i] / ( pow( solution->rho[i], GAMMA) );
		}

		// Region 2
		if ( ( t * u_head ) < (x - x0) && (x - x0) <= t * u_tail ){
			solution->press[i] = pr2;
			solution->rho[i] = rhor2;
			solution->vel[i] = ur2;
			//solution->s[i] = solution->press[i] / ( pow( solution->rho[i], GAMMA) );
		}

		// Region 3
		if ( ( t * u_tail )< ( x - x0 ) && ( x - x0 ) <= t * u_contact ){
			solution->press[i] = pr3;
			solution->rho[i] = rhor3;
			solution->vel[i] = ur3;
			//solution->s[i] = solution->press[i] / ( pow( solution->rho[i], GAMMA) );
		}

		a6 = sqrt(g * pr6/rhor6);

		// Region 4
		if ( t * u_contact < ( x - x0 ) && ( x - x0 ) <= t * u_shock ){
			solution->press[i] = pr4;
			solution->rho[i] = rhor4;
			solution->vel[i] = ur4;
			//solution->s[i] = sol->press[i] / ( pow( solution->rho[i], GAMMA) );
		}

		// Region 5 - Shock!

		// Region 6
		if ( t * u_shock < x - x0 ){
			solution->press[i] = pr6;
			solution->rho[i] = rhor6;
			solution->vel[i] = ur6;
			//solution->s[i] = solution->press[i] / ( pow( solution->rho[i], GAMMA) );
		}
	}

	return;
}

