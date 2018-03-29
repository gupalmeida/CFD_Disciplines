#ifndef INPUT_H
#define INPUT_H

// Output file name
#define filename "output.dat"

// Inputs
#define p4 50.0
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
#define method 3
#define order 1
#define alpha 0.4
#define dissipModel 2

/* Alpha is a constant value used
   in the linear dissipation model.
   The model is selected by the dissipModel
   parameter in which
   0 - linear 2nd order difference model
   1 - linear 4th order difference model
   2 - non-linear jameson model */

/* The order parameter sets the spatial discretization
   order for the vector flux splitting schemes
   which uses one-sided discretization for each
   of the flux vectors (backward for E+ and forward
   for E-).
   The available options are:
         1 - first order
         2 - second order */

/* The solution method is selected by changing
   the method value.
   Available solution methods
   0 - centered scheme
   1 - Lax-Wendroff method
   2 - McCormack method
   3 - Steger and Warming FVS scheme
   4 - 
   5 - 
   6 - 
   7 - 
   8 -  */


#endif
