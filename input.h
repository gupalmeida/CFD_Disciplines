#ifndef INPUT_H
#define INPUT_H

// Output file name
#define filename "output.dat"

// Inputs
#define p4 40.0
#define p1 1.0
#define tol 1.0e-5

// Spatial dicretization parameters
#define l 5.0
#define imax 1001

// Time discretization parameters
#define cfl 0.1
#define tmax 1.0
#define itmax 200000
#define printAt 200

// Initial states
#define U0 0.0
#define T0 288.15

// Fluid constants
#define R 287.0
#define cp 1004.5
#define cv 717.5
#define gamma 1.4

// Method selection
#define method 1

/* Dissipation model parameters */
#define alpha 0.3
#define dissipModel 1
/* For alpha is a constant value used
   in the linear dissipation model.
   The model is selected by the dissipModel
   parameter in which
   0 - linear 2nd order difference model
   1 - linear 4th order difference model
   2 - non-linear jameson model */

#endif
