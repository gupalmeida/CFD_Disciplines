#ifndef INPUT_H
#define INPUT_H

/* Output file name */
#define filename "output.dat"

/* Inputs */
#define p4 90.0
#define p1 1.0
#define tol 1.0e-5

/* Spatial dicretization parameters */
#define l 5.0
#define imax 1001
//#define imax 501

/* Time discretization parameters */
#define cfl 0.1
#define tmax 1.0
#define itmax 200000
#define printAt 200

/* Initial states */
#define U0 0.0
#define T0 288.15

/* Fluid constants */
#define R 287.0
#define cp 1004.5
#define cv 717.5
#define gamma 1.4

/* Method selection */
#define method 7
#define order 2
#define mu 0.4
#define dissipModel 2
#define soundSpeedType 1
#define interfaceSoundSpeed 1
#define limiter 1

#endif

/* The order parameter sets the spatial discretization
   order for the vector flux splitting schemes
   which uses one-sided discretization for each
   of the flux vectors (backward for E+ and forward
   for E-).
   The available options are:
         1 - first order
         2 - second order */

/* mu is a constant value used
   in the linear dissipation model.
   The model is selected by the dissipModel
   parameter in which
   0 - linear 2nd order difference model
   1 - linear 4th order difference model
   2 - non-linear jameson model */

/* soundSpeedType defines the way the code will
   calculate the sound speed at each solution
   point j. The available options are:
       1 - a = aCritical * min( 1.0 , aCritical/|u| )
       2 - a = sqrt( gamma * (p/rho) ) */

/* interfaceSoundSpeed defines the way the code will
   calculate the sound speed at the interface between
   two solution points. The available options are:
       1 - a = min ( a_left, a_right )
       2 - a = 0.5 * ( a_left + a_right )
       3 - a = sqrt( a_left * a_right )
       4 - Roe averaged sound speed */

/* The solution method is selected by changing
   the method value.
   Available solution methods
   0 - centered scheme
   1 - Lax-Wendroff method
   2 - McCormack method
   3 - Steger and Warming FVS scheme
   4 - van Leer FVS non-MUSCL method
   5 - Liou FVS scheme (AUSM+)
   6 - Roe approximate Riemann solver
   7 - Harten TVD method */
