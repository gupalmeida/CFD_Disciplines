#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "input.h"
#include "mesh.h"
#include "aux.h"
#include "solvers.h"

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

    printf("\n\n===============================\n");
    printf("         Centered Scheme\n");
    printf("===============================\n\n");
    
    while (t <= (double) tmax && it < (int) itmax){
        /* calculates the flux vector */
        calcFluxes(solution);

        /* calculates dissipation for the current solution Q */
        calcDissipation(solution);

        /* marching the solution */
        for (int j = 2; j < imax-2; j++){
            Q1np1[j] = solution->Q1[j] - 0.5*(dt/dx)*(solution->E1[j+1] - solution->E1[j-1]) + solution->dissip1[j];
            Q2np1[j] = solution->Q2[j] - 0.5*(dt/dx)*(solution->E2[j+1] - solution->E2[j-1]) + solution->dissip2[j];
            Q3np1[j] = solution->Q3[j] - 0.5*(dt/dx)*(solution->E3[j+1] - solution->E3[j-1]) + solution->dissip3[j];
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
            printf("Iteration: %d, simulation time: %lf, time step: %lf\n",it,t,dt);
        }
    }

    /* deallocating variables */
    free(Q1np1);
    free(Q2np1);
    free(Q3np1);

    printf("\n\n");

}


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

