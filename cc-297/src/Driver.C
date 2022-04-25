#include <iostream>
#include <chrono>
#include <string>

#include "IOobject.H"
#include "BvpOde.H"
//#include "ConvectiveSchemes.H"

double model_prob1_rhs( double x ){ return -1.0; }

int main( int argc, char* argv[] )
{

    SecondOrderODE
    ode_mp1
    (
        1.0,                // coeff of Uxx
        0.0,                // coeff of Ux
        0.0,                // coeff of U
        model_prob1_rhs,    // pointer function to second order ODE
        0.0,                // xMin
        1.0                 // xMax
    );

    BoundaryConditions bc_mp1;
    bc_mp1.setLhsDirichletBc(0.0);
    bc_mp1.setRhsDirichletBc(1.0);

    BvpOde bvpode_mp1( &ode_mp1, &bc_mp1, 6 );
    bvpode_mp1.Solve();

    return 0;
}
