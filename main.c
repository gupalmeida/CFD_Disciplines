#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Prototyping auxiliary functions
double muOperator(double phiNodeJ,double phiNodeJp1);
void stubeAnalytical();
double speedOfSound(double p, double gamma, double rho);


int main(){
    // Method selection
    int method = 1;
    
    // Inputs
    double p4 = 5.0;            // p4/p1
    double p1 = 1.0;

    // Spatial dicretization parameters
    int imax = 1001;                  // # of points - x direction
    double l = 10.0;                  // mesh length
    double dx = (double)(l/(imax-1)); // mesh spacing dx

    // Time discretization parameters
    int itmax = 50;
    double cfl = 1;                   // CFL number
    double tmax = 1.0;                  // maximum simulation time
    double dt = (double) tmax/itmax;    // time step
    
    // Dissipation model
    double alpha = 0.8;
    int dissipModel = 0;

    // ----------------------------------------

    // Initial states
    double U0 = 0.0;            // [m/s]
    double T0 = 288.15;         // [K]

    // Fluid constants
    double R = 287.0;           // [J/kg.K]
    double cp = 1004.5;         // [J/kg.K]
    double cv = 717.5;           // [J/kg.K]
    double gamma = (cp/cv);      // specific heat ratio

    // mesh parameters
    double x[imax];
    
    // equations
    int numberOfEquations = 3;
    double p[imax];
    double q[imax][numberOfEquations],qnp1[imax][numberOfEquations];
    double e[imax][numberOfEquations],enp1[imax][numberOfEquations];
    double eInternal[imax],T[imax],sSpeed[imax];
    double H[imax];
    double dissip[imax][numberOfEquations],jstD[imax][numberOfEquations];
    double J[imax][numberOfEquations];
    
    // output files
    FILE *fp;
     
     
    // =========================================================
    // generating mesh
    x[0] = (double) -l/2.0;
    for (int i = 1; i < imax; i++){
        x[i] = x[i-1] + dx;
        //printf("node: %d, dx: %lf, x: %lf\n",i,dx,x[i]);
    }
    
    // initializing conserved and flux vectors
    double p04 = p4, p01 = p1;
    double rho04 = p4, rho01 = p1;
    double E04 = p04/(gamma-1.0),E01 = p01/(gamma-1.0);
    
    // =========================================================
    //
    //                 INITIALIZING VARIABLES
    //
    // =========================================================
    
    fp = fopen("t0.dat","w+");
    for (int i = 0; i < imax; i++){
        double xPlus = 0.5*(x[i]+fabs(x[i]));
        double xMinus = 0.5*(x[i]-fabs(x[i]));
        p[i] = p04*fabs(xMinus/x[i]) + p01*fabs(xPlus/x[i]);
        // conserved variables vector
        q[i][0] = rho04*fabs(xMinus/x[i]) + rho01*fabs(xPlus/x[i]);
        q[i][1] = 0.0;
        q[i][2] = E04*fabs(xMinus/x[i]) + E01*fabs(xPlus/x[i]);
        qnp1[i][0] = q[i][0];
        qnp1[i][1] = q[i][1];
        qnp1[i][2] = q[i][2];
        // flux vector
        e[i][0] = 0.0;
        e[i][1] = p04*fabs(xMinus/x[i]) + p01*fabs(xPlus/x[i]);
        e[i][2] = 0.0;
        enp1[i][0] = e[i][0];
        enp1[i][1] = e[i][1];
        enp1[i][2] = e[i][2];
        // internal energy
        eInternal[i] = (p[i]/((gamma-1.0)*q[i][0]));
        // temperature
        T[i] = eInternal[i]/cv;
        // speed of sound
        sSpeed[i] = 0.0;
        // jacobian
        J[i][0] = 0.0;
        J[i][1] = 0.0;
        J[i][2] = 0.0;
        // dissipation
        jstD[i][0] = 0.0;
        jstD[i][1] = 0.0;
        jstD[i][2] = 0.0;
        dissip[i][0] = 0.0;
        dissip[i][1] = 0.0;
        dissip[i][2] = 0.0;
        // writting file at t=0 
        fprintf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\n",x[i],q[i][0],(q[i][1]/q[i][0]),eInternal[i],p[i]);
    }
    fclose(fp);
    
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// 
//                      !!!   MARCHING   !!!
// 
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    // time marching parameters
    int it = 0;
    double t = 0.0;

//===================================================================
    // ==========================================================
    //
    //                 BEAM-WARMING METHOD
    //
    // ==========================================================
    //  Q[i]np1 = Q[i]n - 0.5*(dt/dx)*(E[i+1] - E[i-1] + D(Q[i]))
    // ==========================================================
    
    if (method == 1){
        double umax = 0.0;
        while(t != tmax && it != itmax){
            char filename[128] = "";
            for (int i = 0; i < imax; i++){
                double press = (gamma - 1.0)*(q[i][2] - 0.5*(q[i][1]*q[i][1])/q[i][0]);
                double u = q[i][1]/q[i][0];
                if (u > umax){
                    umax = u;
                }
                p[i] = press;
                e[i][0] = q[i][1];
                e[i][1] = q[i][1]*u + press;
                e[i][2] = (q[i][2] + press)*u;
            }
            
            // ARTIFICIAL/NUMERICAL DISSIPATION
            if (dissipModel == 0){
                // 2nd order constant coeff dissipation
/*                
                for (int i = 1; i < imax-1; i++){
                    jstD[i][0] = q[i+1][0] - q[i][0];
                    jstD[i][1] = q[i+1][1] - q[i][1];
                    jstD[i][2] = q[i+1][2] - q[i][2];
                }
                
                for (int i = 1; i < imax-1; i++){
                    dissip[i][0] = alpha*(jstD[i+1][0] - jstD[i-1][0]);
                    dissip[i][1] = alpha*(jstD[i+1][1] - jstD[i-1][1]);
                    dissip[i][2] = alpha*(jstD[i+1][2] - jstD[i-1][2]);
                }
*/
                for (int i = 1; i < imax-1; i++){
                    dissip[i][0] = (alpha/8.0)*(dx/dt)*(q[i+1][0] - 2.0*q[i][0] + q[i-1][0]);
                    dissip[i][1] = (alpha/8.0)*(dx/dt)*(q[i+1][1] - 2.0*q[i][1] + q[i-1][1]);
                    dissip[i][2] = (alpha/8.0)*(dx/dt)*(q[i+1][2] - 2.0*q[i][2] + q[i-1][2]);
                }
            }
            else if (dissipModel == 1){
                // 4th order constant coeff dissipation
/*
                for (int i = 1; i < imax-2; i++){
                    jstD[i][0] = q[i+2][0] - 3.0*q[i+1][0] + 3.0*q[i][0] - q[i-1][0];
                    jstD[i][1] = q[i+2][1] - 3.0*q[i+1][1] + 3.0*q[i][1] - q[i-1][1];
                    jstD[i][2] = q[i+2][2] - 3.0*q[i+1][2] + 3.0*q[i][2] - q[i-1][2];
                }
                
                for (int i = 1; i < imax-1; i++){
                    dissip[i][0] = alpha*(jstD[i-1][0] - jstD[i+1][0]);
                    dissip[i][1] = alpha*(jstD[i-1][1] - jstD[i+1][1]);
                    dissip[i][2] = alpha*(jstD[i-1][2] - jstD[i+1][2]);
                }
*/
                for (int i = 1; i < imax-1; i++){
                    dissip[i][0] = -(alpha/8.0)*(dx/dt)*(q[i+2][0] - 4.0*q[i+1][0] + 6.0*q[i][0] - 4.0*q[i-1][0] + q[i-2][0]);
                    dissip[i][1] = -(alpha/8.0)*(dx/dt)*(q[i+2][1] - 4.0*q[i+1][1] + 6.0*q[i][1] - 4.0*q[i-1][1] + q[i-2][1]);
                    dissip[i][2] = -(alpha/8.0)*(dx/dt)*(q[i+2][2] - 4.0*q[i+1][2] + 6.0*q[i][2] - 4.0*q[i-1][2] + q[i-2][2]);
                }
            }
            dt = cfl*dx/umax;
            // =========================================
        
            for (int i = 1; i < imax-1; i++){
                qnp1[i][0] = q[i][0] - 0.5*(dt/dx)*(e[i+1][0] - e[i-1][0] + dissip[i][0]);
                qnp1[i][1] = q[i][1] - 0.5*(dt/dx)*(e[i+1][1] - e[i-1][1] + dissip[i][1]);
                qnp1[i][2] = q[i][2] - 0.5*(dt/dx)*(e[i+1][2] - e[i-1][2] + dissip[i][1]);
            }

            for (int i = 0; i < imax; i++){
                q[i][0] = qnp1[i][0];
                q[i][1] = qnp1[i][1];
                q[i][2] = qnp1[i][2];
                eInternal[i] = (p[i]/((gamma-1.0)*q[i][0]));
                T[i] = eInternal[i]/cv;
            }
            
            if (it%2 == 0){
                sprintf(filename,"iter%d.dat",it);
                fp = fopen(filename,"w+");
                for (int i = 0; i < imax; i++){
                    fprintf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\n",x[i],q[i][0],(q[i][1]/q[i][0]),eInternal[i],p[i]);
                }
                fclose(fp);
            }
            t = t + dt;
            it++;
       }
    }
//===================================================================
    // =========================================================
    //
    //                 LAX-WENDROFF METHOD
    //
    // =========================================================
    //  Q[i]np1 = Q[i]n - 0.5*(dt/dx)*(E[i+1] - E[i-1]) //
    //          + 0.5*(dt²/dx²)(J[i+0.5]*(E[i+1] - E[i])//
    //          - J[i-0.5]*(E[i] - E[i-1]))
    // =========================================================
    
    if (method == 2){
        double diff1[imax][3], diff2[imax][3];
        // initializing auxiliary vectors
        for (int i = 0; i < imax; i++){
            diff1[i][0] = 0.0;
            diff1[i][1] = 0.0;
            diff1[i][2] = 0.0;
            diff2[i][0] = 0.0;
            diff2[i][1] = 0.0;
            diff2[i][2] = 0.0;
            printf("dissip: %lf, %lf, %lf\n",dissip[i][0],dissip[i][2],dissip[i][2]);
        }
        
        //while(t != tmax && it != itmax){
        while(it != 200){
            t = t + dt;
            it++;
            char filename[128] = "";
            for (int i = 0; i < imax; i++){
                p[i] = (gamma - 1.0)*(q[i][2] - 0.5*pow(q[i][1],2.0)/q[i][0]);
                e[i][0] = q[i][1];
                e[i][1] = q[i][1]*(q[i][1]/q[i][0]) + p[i];
                e[i][2] = (q[i][2] + p[i])*(q[i][1]/q[i][0]);
            }
            
            for (int i = 1; i < imax-1; i++){
                // (E[j+1] - E[j])
                diff1[i][0] = e[i+1][0] - e[i][0];
                diff1[i][1] = e[i+1][1] - e[i][1];
                diff1[i][2] = e[i+1][2] - e[i][2];

                // (E[j] - E[j-1])
                diff2[i][0] = e[i][0] - e[i-1][0];
                diff2[i][1] = e[i][1] - e[i-1][1];
                diff2[i][2] = e[i][2] - e[i-1][2];
            }

            for (int i = 1; i < imax-1; i++){
                // common for both j+half and j-half
                double ujm1 = (q[i-1][1]/q[i-1][0]);
                double uj = (q[i][1]/q[i][0]);
                double ujp1 = (q[i+1][1]/q[i+1][0]);
                double rhojm1 = q[i-1][0];
                double rhoj = q[i][0];
                double rhojp1 = q[i+1][0];
                double ejm1 = q[i-1][2];
                double ej = q[i][2];
                double ejp1 = q[i+1][2];
                // j + half
                // ==================
                double a11half = 0.0;
                double a12half = 1.0;
                double a13half = 0.0;
                // ==================
                double a21phalf = 0.5*((gamma-3.0)/2.0)*(pow(uj,2.0) + pow(ujp1,2.0));
                double a22phalf = 0.5*(3.0-gamma)*(uj + ujp1);
                double a23half = (gamma-1.0);
                // ==================
                double a31phalf = 0.5*(((gamma-1.0)*pow(uj,3.0) - gamma*((uj/rhoj)*ej)) + ((gamma-1.0)*pow(ujp1,3.0) - gamma*((ujp1/rhojp1)*ejp1)));
                double a32phalf = 0.5*(((gamma/rhoj)*ej - 1.5*(gamma-1.0)*(pow(uj,3.0))) + ((gamma/rhojp1)*ejp1 - 1.5*(gamma-1.0)*(pow(ujp1,3.0))));
                double a33phalf = 0.5*gamma*(uj + ujp1);
                // ==================
                // j - half
                double a21mhalf = 0.5*((gamma-3.0)/2.0)*(pow(uj,2.0) + pow(ujm1,2.0));
                double a22mhalf = 0.5*(3.0-gamma)*(uj + ujm1);
                // ==================
                double a31mhalf = 0.5*(((gamma-1.0)*pow(uj,3.0) - gamma*((uj/rhoj)*ej)) + ((gamma-1.0)*pow(ujm1,3.0) - gamma*((ujm1/rhojm1)*ejm1)));
                double a32mhalf = 0.5*(((gamma/rhoj)*ej - 1.5*(gamma-1.0)*(pow(uj,3.0))) + ((gamma/rhojp1)*ejm1 - 1.5*(gamma-1.0)*(pow(ujm1,3.0))));
                double a33mhalf = 0.5*gamma*(uj + ujm1);
                // --------------------------------------------------------------------
                // --------------------------------------------------------------------
                // --------------------------------------------------------------------
                dissip[i][0] = (0.5*pow((dt/dx),2.0))*(a11half*diff1[i][0]+a12half*diff1[i][1]+a13half*diff1[i][2] - (a11half*diff2[i][0]+a12half*diff2[i][1]+a13half*diff2[i][2]));
                dissip[i][1] = (0.5*pow((dt/dx),2.0))*(a21phalf*diff1[i][0]+a22phalf*diff1[i][1]+a23half*diff1[i][2] - (a21mhalf*diff2[i][0]+a22mhalf*diff2[i][1]+a23half*diff2[i][2]));
                dissip[i][2] = (0.5*pow((dt/dx),2.0))*(a31phalf*diff1[i][0]+a32phalf*diff1[i][1]+a33phalf*diff1[i][2] - (a31mhalf*diff2[i][0]+a32mhalf*diff2[i][1]+a33mhalf*diff2[i][2]));
                printf("x: %lf, d[0]: %lf, d[1]: %lf, d[2]: %lf\n",x[i],diff1[i][0],diff1[i][1],diff1[i][2]);
            }

    
            for (int i = 1; i < imax-1; i++){
                double lbd1 = 0.5*(dt/dx);
                double lbd2 = 0.5*(dt*dt)/(dx*dx);
                qnp1[i][0] = q[i][0] - lbd1*(e[i+1][0] - e[i-1][0]) + lbd2*dissip[i][0];
                qnp1[i][1] = q[i][1] - lbd1*(e[i+1][1] - e[i-1][1]) + lbd2*dissip[i][1];
                qnp1[i][2] = q[i][2] - lbd1*(e[i+1][2] - e[i-1][2]) + lbd2*dissip[i][2];
            }
            
            for (int i = 0; i < imax; i++){
                q[i][0] = qnp1[i][0];
                q[i][1] = qnp1[i][1];
                q[i][2] = qnp1[i][2];
                eInternal[i] = (p[i]/((gamma-1.0)*q[i][0]));
                T[i] = eInternal[i]/cv;
            }
            
            if (it%50 == 0){
                sprintf(filename,"iter%d.dat",it);
                fp = fopen(filename,"w+");
                for (int i = 0; i < imax; i++){
                    fprintf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\n",x[i],q[i][0],(q[i][1]/q[i][0]),eInternal[i],p[i]);
                }
                fclose(fp);
            }
        }
    }
//===================================================================
    // =========================================================
    //
    //                 MAcCORMACK METHOD
    //
    // =========================================================
//===================================================================
    // =========================================================
    //
    //                 ---------- METHOD
    //
    // =========================================================
//===================================================================
    // =========================================================
    //
    //                 ---------- METHOD
    //
    // =========================================================
//===================================================================
    // =========================================================
    //
    //                 ---------- METHOD
    //
    // =========================================================
//===================================================================
    // =========================================================
    //
    //                 ---------- METHOD
    //
    // =========================================================
//===================================================================
    // =========================================================
    //
    //                 ---------- METHOD
    //
    // =========================================================
//===================================================================
    // =========================================================
    //
    //                 ---------- METHOD
    //
    // =========================================================









}
//===================================================================
//===================================================================
//===================================================================
//===================================================================
//===================================================================
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// 
//                      AUXILIARY FUNCTIONS
// 
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

double muOperator(double phiNodeJ,double phiNodeJp1){
    // This function calculates the mean value of a given property at the interface
    double mu;
    mu = 0.5*(phiNodeJ + phiNodeJp1);
    return mu;
}

double speedOfSound(double p, double gamma, double rho){
    double a;
    a = sqrt(gamma*p/rho);
    return a;
}

void analytical(){
    double p1,p2,p3,p4,p5;
    double r1,r2,r3,r4,r5;
    double u1,u2,u3,u4,u5;
    double gm = 1.4;
    double c1,c2,c3;
    double a5,b5;
    double tmax = 1.0;
    double x,x0 = 0.0;
    double uTail,uHead,uContact,uShock;
    double f,flinha,p0,pOld,pNew;
}


