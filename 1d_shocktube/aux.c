#include "aux.h"

void allocSolution(results * solution){
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
        solution->E1p[i] = 0.0;
        solution->E2p[i] = 0.0;
        solution->E3p[i] = 0.0;
        solution->E1m[i] = 0.0;
        solution->E2m[i] = 0.0;
        solution->E3m[i] = 0.0;
    }

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

void calcVanLeerFluxes(results * solution){
    double a, u, p, rho, M;
    double fmassp, fmassm;
    double gm = gamma;

    for (int j = 0; j < imax; j++){
        rho = solution->rho[j];
        u = solution->vel[j];
        p = solution->press[j];

        a = sqrt(gamma*(p/rho));

        M = u/a;
        fmassp = (rho*a/4.0)*(M + 1.0)*(M + 1.0);
        fmassm = -(rho*a/4.0)*(M - 1.0)*(M - 1.0);

        if (M <= -1.0){
            solution->E1p[j] = 0.0;
            solution->E2p[j] = 0.0;
            solution->E3p[j] = 0.0;

            solution->E1m[j] = solution->Q2[j];
            solution->E2m[j] = (solution->Q2[j] * solution->Q2[j])/solution->Q1[j] + solution->press[j];
            solution->E3m[j] = (solution->Q3[j] + solution->press[j])*solution->vel[j];
        }
        else if(M >= 1.0){
            solution->E1m[j] = 0.0;
            solution->E2m[j] = 0.0;
            solution->E3m[j] = 0.0;

            solution->E1p[j] = solution->Q2[j];
            solution->E2p[j] = (solution->Q2[j] * solution->Q2[j])/solution->Q1[j] + solution->press[j];
            solution->E3p[j] = (solution->Q3[j] + solution->press[j])*solution->vel[j];
        }
        else {
            solution->E1p[j] = fmassp;
            solution->E2p[j] = fmassp*(1.0/gm)*((gm-1.0)*u + 2.0*a);
            solution->E3p[j] = fmassp*(
            ((gm-1.0)*u + 2.0*a)*((gm-1.0)*u + 2.0*a)
            /(2.0*(gm*gm-1.0)));

            solution->E1m[j] = fmassm;
            solution->E2m[j] = fmassm*(1.0/gm)*((gm-1.0)*u - 2.0*a);
            solution->E3m[j] = fmassm*(
            ((gm-1.0)*u - 2.0*a)*((gm-1.0)*u - 2.0*a)
            /(2.0*(gm*gm-1.0)));
        }
    }
    //printf("l1: %f, l2: %f, l3: %f\n",lbd1,lbd2,lbd3);

    return;
}

void calcLiouFluxes( results * solution ){

    double phiL[3], phiR[3];
    double hL, hR;
    double rhoL, rhoR;

    double aL, aR;
	double aCriticalL, aCriticalR;
    double ajph;

    double uL, uR;
    double pjph;
    double pL, pLplus, pR, pRminus;
    double pAlphaPlus, pAlphaMinus;

    double ML, MR;
    double MLplus, MRminus;
    double MbetaPlus, MbetaMinus;
    double mach;
    //double machPlus, machMinus;

    double alpha = 3.0/16.0;
    double beta = 1.0/8.0;

    for (int j = 1; j < imax-2; j++){
        uL = solution->vel[j];
        uR = solution->vel[j+1];
        pL = solution->press[j];
        pR = solution->press[j+1];
        rhoL = solution->Q1[j];
        rhoR = solution->Q1[j+1];
        hL = (solution->Q3[j]+solution->press[j])/solution->Q1[j];
        hR = (solution->Q3[j+1]+solution->press[j+1])/solution->Q1[j+1];

        /* computing phiL and phiR */
        phiL[0] = solution->Q1[j];
        phiL[1] = solution->Q2[j];
        phiL[2] = solution->Q3[j] + pL;

        phiR[0] = solution->Q1[j+1];
        phiR[1] = solution->Q2[j+1];
        phiR[2] = solution->Q3[j+1] + pR;

        /* computing the sound speed at interface j+1/2 */
        aCriticalL = sqrt( (2.0*(gamma-1.0)/(gamma+1.0))*hL );
        aCriticalR = sqrt( (2.0*(gamma-1.0)/(gamma+1.0))*hR );

        switch (soundSpeedType){
            case 1:
                aL = aCriticalL * min( 1.0, (aCriticalL/fabs(uL)) );
                aR = aCriticalR * min( 1.0, (aCriticalR/fabs(uR)) );
                break;
            case 2:
                aL = sqrt( gamma* (pL/rhoL) );
                aR = sqrt( gamma* (pR/rhoR) );
                break;
            default:
                aL = aCriticalL * min( 1.0, (aCriticalL/fabs(uL)) );
                aR = aCriticalR * min( 1.0, (aCriticalR/fabs(uR)) );
                break;
        }

        switch (interfaceSoundSpeed){
            case 1:
                ajph = min( aL,aR );
                break;
            case 2:
                ajph = 0.5*( aL + aR );
                break;
            case 3:
                ajph = sqrt( aL * aR );
                break;
            case 4:
                //reserved for roe averages
                break;
            default:
                ajph = min( aL,aR );
                break;
        }

        /* computing (M,p)j+1/2 */
        ML = fabs( uL ) / ajph;
        MR = fabs( uR ) / ajph;

        MbetaPlus = 0.25*(ML + 1.0)*(ML + 1.0)
               + beta*(ML*ML - 1.0)*(ML*ML - 1.0);
        MbetaMinus = -0.25*(MR - 1.0)*(MR - 1.0)
               - beta*(MR*MR - 1.0)*(MR*MR - 1.0);

        pAlphaPlus = 0.25*(ML + 1.0)*(ML + 1.0)*(2.0 - ML)
               + alpha*ML*(ML*ML - 1.0)*(ML*ML - 1.0);
        pAlphaMinus = 0.25*(MR - 1.0)*(MR - 1.0)*(2.0 + MR)
               - alpha*MR*(MR*MR - 1.0)*(MR*MR - 1.0);

        if ( fabs(ML) >= 1.0 ){
            MLplus = 0.5*( ML + fabs(ML) );
            pLplus = 0.5*( 1.0 + (ML / fabs(ML)) );
        }
        else{
            MLplus = MbetaPlus;
            pLplus = pAlphaPlus;
        }

        if ( fabs(MR) >= 1.0 ){
            MRminus = 0.5*( MR - fabs(MR) );
            pRminus = 0.5*( 1.0 - (MR / fabs(MR)) );
        }
        else{
            MRminus = MbetaMinus;
            pRminus = pAlphaMinus;
        }

        mach = MLplus + MRminus;
        //machPlus = 0.5*( mach + fabs(mach) );
        //machMinus = 0.5*( mach - fabs(mach) );
        pjph = pLplus*pL + pRminus*pR;
        
        /* computing the fluxes at j+1/2 */
        solution->E1p[j] = 0.5 * ajph * mach * (phiL[0] + phiR[0])
           - 0.5 * ajph * fabs( mach )*( phiR[0] - phiL[0] );
        solution->E2p[j] = 0.5 * ajph * mach * (phiL[1] + phiR[1])
           - 0.5 * ajph * fabs( mach )*( phiR[1] - phiL[1] ) + pjph;
        solution->E3p[j] = 0.5 * ajph * mach * (phiL[2] + phiR[2])
           - 0.5 * ajph * fabs( mach )*( phiR[2] - phiL[2] );

    }

	return;
}

void calcRoeFluxes( results * solution ){

    double uL, uR;
    double pL, pR;
    double rhoL, rhoR;
    double hL, hR;
    double uRoe, rhoRoe, hRoe, aRoe;

    double lbd1, lbd2, lbd3;
    double r1[3], r2[3], r3[3];
    double deltaQ[3];
    double dp, du, drho;
    double alpha1, alpha2, alpha3;

    for (int j = 1; j < imax-2; j++){
        uL = solution->vel[j];
        uR = solution->vel[j+1];
        pL = solution->press[j];
        pR = solution->press[j+1];
        rhoL = solution->rho[j];
        rhoR = solution->rho[j+1];
        hL = (solution->Q3[j] + solution->press[j])/solution->Q1[j];
        hR = (solution->Q3[j+1] + solution->press[j+1])/solution->Q1[j+1];
        du = uR - uL;
        dp = pR - pL;
        drho = rhoR - rhoL;

        /* computing roe averaged variables */
        rhoRoe = sqrt( rhoL * rhoR );
        uRoe = (sqrt( rhoL ) * uL + sqrt( rhoR ) * uR)
               / ( sqrt( rhoL ) + sqrt( rhoR ) );
        hRoe = (sqrt( rhoL ) * hL + sqrt( rhoR ) * hR)
               / ( sqrt( rhoL ) + sqrt( rhoR ) );
        aRoe = sqrt( (gamma-1.0)*(hRoe - 0.5 * uRoe * uRoe) );

        /* computing the lambda eigenvalues */
        lbd1 = uRoe - aRoe;
        lbd2 = uRoe;
        lbd3 = uRoe + aRoe;

        /* computing the rk eigenvectors */
        r1[0] = 1.0;
        r1[1] = lbd1;
        r1[2] = hRoe - uRoe*aRoe;
        r2[0] = 1.0;
        r2[1] = lbd2;
        r2[2] = 0.5*uRoe*uRoe;
        r3[0] = 1.0;
        r3[1] = lbd3;
        r3[2] = hRoe + uRoe*aRoe;

        /* computing characteristic properties */
        //alpha1 = du - (1.0/(rhoRoe * aRoe))*dp;
        //alpha2 = drho - (1.0/(aRoe * aRoe))*dp;
        //alpha3 = du + (1.0/(rhoRoe * aRoe))*dp;
        alpha1 = 0.5 * (dp - rhoRoe * aRoe * du) / (aRoe*aRoe);
        alpha2 = - (dp - drho * aRoe * aRoe) / (aRoe*aRoe);
        alpha3 = 0.5 * (dp + rhoRoe * aRoe * du) / (aRoe*aRoe);

        /* computing step in Q */
        deltaQ[0] = fabs( lbd1 ) * alpha1 * r1[0]
                  + fabs( lbd2 ) * alpha2 * r2[0]
                  + fabs( lbd3 ) * alpha3 * r3[0];
        deltaQ[1] = fabs( lbd1 ) * alpha1 * r1[1]
                  + fabs( lbd2 ) * alpha2 * r2[1]
                  + fabs( lbd3 ) * alpha3 * r3[1];
        deltaQ[2] = fabs( lbd1 ) * alpha1 * r1[2]
                  + fabs( lbd2 ) * alpha2 * r2[2]
                  + fabs( lbd3 ) * alpha3 * r3[2];

        /* computing the fluxes */
        solution->E1p[j] = 0.5 * ( solution->E1[j] + solution->E1[j+1] - deltaQ[0]);
        solution->E2p[j] = 0.5 * ( solution->E2[j] + solution->E2[j+1] - deltaQ[1]);
        solution->E3p[j] = 0.5 * ( solution->E3[j] + solution->E3[j+1] - deltaQ[2]);

    }

}

void calcHartenFluxes( results * solution, double lambda ){

    double uL, uR;
    //double pL, pR;
    double rhoL, rhoR;
    double eL, eR;
    double hL, hR;
    double uRoe, hRoe, aRoe;

    double lbd1, lbd2, lbd3;
    double r1[3], r2[3], r3[3];
    double deltaQ[3];
    //double dp, du;
    double drho, de, dqm;
    double alpha1, alpha2, alpha3;
    double C1, C2;
    double g[3][imax], gBar[3][imax];
    double nu[3], ksi[3], s[3], eps[3], gm[3];

    eps[0] = 0.05;
    eps[1] = 0.0;
    eps[2] = 0.0;

    for (int j = 1; j < imax-2; j++){
        uL = solution->vel[j];
        uR = solution->vel[j+1];
        //pL = solution->press[j];
        //pR = solution->press[j+1];
        rhoL = solution->rho[j];
        rhoR = solution->rho[j+1];
        eL = solution->Q3[j];
        eR = solution->Q3[j+1];
        hL = (solution->Q3[j] + solution->press[j])/solution->Q1[j];
        hR = (solution->Q3[j+1] + solution->press[j+1])/solution->Q1[j+1];
        //du = uR - uL;
        //dp = pR - pL;
        drho = rhoR - rhoL;
        de = eR - eL;
        dqm = solution->Q2[j+1] - solution->Q2[j];

        /* computing roe averaged variables */
        uRoe = (sqrt( rhoL ) * uL + sqrt( rhoR ) * uR)
               / ( sqrt( rhoL ) + sqrt( rhoR ) );
        hRoe = (sqrt( rhoL ) * hL + sqrt( rhoR ) * hR)
               / ( sqrt( rhoL ) + sqrt( rhoR ) );
        aRoe = sqrt( (gamma-1.0)*(hRoe - 0.5 * uRoe * uRoe) );

        /* computing the lambda eigenvalues */
        lbd1 = uRoe - aRoe;
        lbd2 = uRoe;
        lbd3 = uRoe + aRoe;

        /* computing characteristic properties */
        C1 = ((gamma-1.0)/(aRoe*aRoe))*(de
                + 0.5 * uRoe * uRoe * drho
                - uRoe * dqm);
        C2 = (1.0/aRoe) * (dqm - uRoe * drho);

        alpha1 = 0.5 * ( C1 - C2 );
        alpha2 = drho - C1;
        alpha3 = 0.5 * ( C1 + C2 );

        /* computing step in Q */
        
        nu[0] = lambda * lbd1;
        nu[1] = lambda * lbd2;
        nu[2] = lambda * lbd3;

        switch (order){
            case 1:
                gBar[0][j] = 0.0;
                gBar[1][j] = 0.0;
                gBar[2][j] = 0.0;
                break;
            case 2:
                ksi[0] = ksiHarten(nu[0],eps[0]);
                ksi[1] = ksiHarten(nu[1],eps[1]);
                ksi[2] = ksiHarten(nu[2],eps[2]);

                gBar[0][j] = 0.5 * ( ksi[0] - nu[0]*nu[0] ) * alpha1;
                gBar[1][j] = 0.5 * ( ksi[1] - nu[1]*nu[1] ) * alpha2;
                gBar[2][j] = 0.5 * ( ksi[2] - nu[2]*nu[2] ) * alpha3;
                break;
            default:
                gBar[0][j] = 0.0;
                gBar[1][j] = 0.0;
                gBar[2][j] = 0.0;
                break;
        }
    }

    for (int j = 1; j < imax-2; j++){
        switch (order){
            case 1:
                g[0][j] = 0.0;
                g[1][j] = 0.0;
                g[2][j] = 0.0;
                break;
            case 2:
                switch (limiter){
                    case 0:
                        s[0] = copysign( 1.0, gBar[0][j] );
                        s[1] = copysign( 1.0, gBar[1][j] );
                        s[2] = copysign( 1.0, gBar[2][j] );

                        g[0][j] = gBar[0][j];
                        g[1][j] = gBar[1][j];
                        g[2][j] = gBar[2][j];
                        break;
                    case 1:
                        s[0] = copysign( 1.0, gBar[0][j] );
                        s[1] = copysign( 1.0, gBar[1][j] );
                        s[2] = copysign( 1.0, gBar[2][j] );

                        g[0][j] =s[0]*max( 0.0, min(fabs(gBar[0][j]),s[0]*gBar[0][j-1]) );
                        g[1][j] =s[1]*max( 0.0, min(fabs(gBar[1][j]),s[1]*gBar[1][j-1]) );
                        g[2][j] =s[2]*max( 0.0, min(fabs(gBar[2][j]),s[2]*gBar[2][j-1]) );
                        break;
                    default:
                        s[0] = copysign( 1.0, gBar[0][j] );
                        s[1] = copysign( 1.0, gBar[1][j] );
                        s[2] = copysign( 1.0, gBar[2][j] );

                        g[0][j] =s[0]*max( 0.0, min(fabs(gBar[0][j]),s[0]*gBar[0][j-1]) );
                        g[1][j] =s[1]*max( 0.0, min(fabs(gBar[1][j]),s[1]*gBar[1][j-1]) );
                        g[2][j] =s[2]*max( 0.0, min(fabs(gBar[2][j]),s[2]*gBar[2][j-1]) );
                        break;
                }
                break;
            default:
                g[0][j] = 0.0;
                g[1][j] = 0.0;
                g[2][j] = 0.0;
                break;
        }
    }

    for (int j = 1; j < imax-2; j++){
        uL = solution->vel[j];
        uR = solution->vel[j+1];
        //pL = solution->press[j];
        //pR = solution->press[j+1];
        rhoL = solution->rho[j];
        rhoR = solution->rho[j+1];
        eL = solution->Q3[j];
        eR = solution->Q3[j+1];
        hL = (solution->Q3[j] + solution->press[j])/solution->Q1[j];
        hR = (solution->Q3[j+1] + solution->press[j+1])/solution->Q1[j+1];
        //du = uR - uL;
        //dp = pR - pL;
        drho = rhoR - rhoL;
        de = eR - eL;
        dqm = solution->Q2[j+1] - solution->Q2[j];

        /* computing roe averaged variables */
        //rhoRoe = sqrt( rhoL * rhoR );
        uRoe = (sqrt( rhoL ) * uL + sqrt( rhoR ) * uR)
               / ( sqrt( rhoL ) + sqrt( rhoR ) );
        hRoe = (sqrt( rhoL ) * hL + sqrt( rhoR ) * hR)
               / ( sqrt( rhoL ) + sqrt( rhoR ) );
        aRoe = sqrt( (gamma-1.0)*(hRoe - 0.5 * uRoe * uRoe) );

        /* computing the lambda eigenvalues */
        lbd1 = uRoe - aRoe;
        lbd2 = uRoe;
        lbd3 = uRoe + aRoe;

        /* computing the rk eigenvectors */
        r1[0] = 1.0;
        r1[1] = lbd1;
        r1[2] = hRoe - uRoe*aRoe;
        r2[0] = 1.0;
        r2[1] = lbd2;
        r2[2] = 0.5*uRoe*uRoe;
        r3[0] = 1.0;
        r3[1] = lbd3;
        r3[2] = hRoe + uRoe*aRoe;

        /* computing characteristic properties */
        C1 = ((gamma-1.0)/(aRoe*aRoe))*(de
                + 0.5 * uRoe * uRoe * drho
                - uRoe * dqm);
        C2 = (1.0/aRoe) * (dqm - uRoe * drho);

        alpha1 = 0.5 * ( C1 - C2 );
        alpha2 = drho - C1;
        alpha3 = 0.5 * ( C1 + C2 );

        /* computing gamma k */
        if (fabs( alpha1 ) > eps[0]){
            gm[0] = (g[0][j+1] - g[0][j]) / alpha1;
        }
        else {
            gm[0] = 0.0;
        }
        if (fabs( alpha2 ) > eps[1]){
            gm[1] = (g[1][j+1] - g[1][j]) / alpha2;
        }
        else {
            gm[1] = 0.0;
        }
        if (fabs( alpha3 ) > eps[2]){
            gm[2] = (g[2][j+1] - g[2][j]) / alpha3;
        }
        else {
            gm[2] = 0.0;
        }

        /* computing step in Q */
        nu[0] = lambda * lbd1;
        nu[1] = lambda * lbd2;
        nu[2] = lambda * lbd3;

        ksi[0] = ksiHarten( (nu[0] + gm[0]) , eps[0] );
        ksi[1] = ksiHarten( (nu[1] + gm[1]) , eps[1] );
        ksi[2] = ksiHarten( (nu[2] + gm[2]) , eps[2] );

        deltaQ[0] = r1[0] * (g[0][j] + g[0][j+1] - ksi[0] * alpha1)
                  + r2[0] * (g[1][j] + g[1][j+1] - ksi[1] * alpha2)
                  + r3[0] * (g[2][j] + g[2][j+1] - ksi[2] * alpha3);
        deltaQ[1] = r1[1] * (g[0][j] + g[0][j+1] - ksi[0] * alpha1)
                  + r2[1] * (g[1][j] + g[1][j+1] - ksi[1] * alpha2)
                  + r3[1] * (g[2][j] + g[2][j+1] - ksi[2] * alpha3);
        deltaQ[2] = r1[2] * (g[0][j] + g[0][j+1] - ksi[0] * alpha1)
                  + r2[2] * (g[1][j] + g[1][j+1] - ksi[1] * alpha2)
                  + r3[2] * (g[2][j] + g[2][j+1] - ksi[2] * alpha3);

        /* computing the fluxes */
        solution->E1p[j] = 0.5 * ( solution->E1[j] + solution->E1[j+1] + (1.0/lambda)*deltaQ[0]);
        solution->E2p[j] = 0.5 * ( solution->E2[j] + solution->E2[j+1] + (1.0/lambda)*deltaQ[1]);
        solution->E3p[j] = 0.5 * ( solution->E3[j] + solution->E3[j+1] + (1.0/lambda)*deltaQ[2]);

    }

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
                solution->dissip1[j] = (mu/8.0)*(1.0/lambda)*(solution->Q1[j+1] - 2.0*solution->Q1[j] + solution->Q1[j-1]);
                solution->dissip2[j] = (mu/8.0)*(1.0/lambda)*(solution->Q2[j+1] - 2.0*solution->Q2[j] + solution->Q2[j-1]);
                solution->dissip3[j] = (mu/8.0)*(1.0/lambda)*(solution->Q3[j+1] - 2.0*solution->Q3[j] + solution->Q3[j-1]);
            }
            break;

        case 1:

            for (int j = 2; j < imax-2; j++){
                solution->dissip1[j] = -(mu/8.0)*(1.0/lambda)*(solution->Q1[j+2] - 4.0*solution->Q1[j+1] + 6.0*solution->Q1[j] - 4.0*solution->Q1[j-1] + solution->Q1[j-2]);
                solution->dissip2[j] = -(mu/8.0)*(1.0/lambda)*(solution->Q2[j+2] - 4.0*solution->Q2[j+1] + 6.0*solution->Q2[j] - 4.0*solution->Q2[j-1] + solution->Q2[j-2]);
                solution->dissip3[j] = -(mu/8.0)*(1.0/lambda)*(solution->Q3[j+2] - 4.0*solution->Q3[j+1] + 6.0*solution->Q3[j] - 4.0*solution->Q3[j-1] + solution->Q3[j-2]);
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
                solution->dissip1[j] = (mu/8.0)*(1.0/lambda)*(solution->Q1[j+1] - 2.0*solution->Q1[j] + solution->Q1[j-1]);
                solution->dissip2[j] = (mu/8.0)*(1.0/lambda)*(solution->Q2[j+1] - 2.0*solution->Q2[j] + solution->Q2[j-1]);
                solution->dissip3[j] = (mu/8.0)*(1.0/lambda)*(solution->Q3[j+1] - 2.0*solution->Q3[j] + solution->Q3[j-1]);
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

double min(double a, double b){
    if (a < b){
        return a;
    }
    else{
        return b;
    }
}

double ksiHarten( double a, double eps ){
    double ksi;

    if ( fabs( a ) < eps ){
        ksi = 0.5 * ( ((a*a)/eps) + eps );
    }
    else {
        ksi = fabs( a );
    }

    return ksi;
}
