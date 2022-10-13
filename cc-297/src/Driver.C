#include <iostream>
#include <chrono>
#include <string>

#include "IOobject.H"
#include "BvpOde.H"
#include "Aux.H"
//#include "ConvectiveSchemes.H"

double model_prob1_rhs( double x ){ return -std::sin(x); }

int main( int argc, char* argv[] )
{
    int imax = readInput<int>( "IMAX", SETUP_FILE );

    SecondOrderODE
    ode_mp1
    (
        1.0,                // coeff of Uxx
        0.0,                // coeff of Ux
        0.0,                // coeff of U
        model_prob1_rhs,    // pointer function to second order ODE
        0.0,                // xMin
        PI                  // xMax
    );

    BoundaryConditions bc_mp1;
    bc_mp1.setLhsDirichletBc(0.0);
    bc_mp1.setRhsDirichletBc(1.0);

    BvpOde bvpode_mp1( &ode_mp1, &bc_mp1, imax );
    bvpode_mp1.Solve();

    return 0;
}
