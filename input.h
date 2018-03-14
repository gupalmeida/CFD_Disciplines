#ifndef INPUT_H
#define INPUT_H

// Output file name
#define filename "output.dat"

// Inputs
#define p4 5.0
#define p1 1.0

// Spatial dicretization parameters
//double l = 10.0;
//int imax = 1001;
#define l 10.0
#define imax 21
//double dx (l/(imax-1)) // mesh spacing dx


// Time discretization parameters
//double cfl = 0.5;
//double tmax = 1.0;
//double dt = 5.0e-3;
//double itmax = 200;
#define cfl 0.5
#define tmax 1.0
#define itmax 200

// Initial states
//double U0 = 0.0;
//double T0 = 288.15;
#define U0 0.0
#define T0 288.15

// Fluid constants
//double R = 287.0;
//double cp = 1004.5;
//double cv = 717.5;
#define R 287.0
#define cp 1004.5
#define cv 717.5
#define gamma 1.4

// Method selection
//int method = 1;
#define method 1

// Dissipation model
//double alpha = 0.8;
//int dissipModel = 0;
#define alpha 0.8
#define dissipModel 0

#endif
