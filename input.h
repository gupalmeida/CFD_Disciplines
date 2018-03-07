// Inputs
#define p4 5            // p4/p1 = 5
#define p1 1

// Spatial dicretization parameters
#define l 10            // mesh length
#define imax 1001       // number of mesh points - x direction
#define dx (l/(imax-1)) // mesh spacing dx

// Time discretization parameters
#define cfl 0.5         // Courant-Friedrichs-Lewy number
#define tmax 1.0        // maximum simulation time
#define dt 5.0e-3       // time step
#define itmax 200

// ----------------------------------------

// Initial states
#define U0 0            // [m/s]
#define T0 288.15       // [K]

// Fluid constants
#define R 287           // [J/kg.K]
#define cp 1004.5       // [J/kg.K]
#define cv 717.5        // [J/kg.K]
#define gama cp/cv      // specific heat ratio

// Method selection
#define method 1
