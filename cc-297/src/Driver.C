#include <iostream>
#include <chrono>
#include <string>

#include "IOobject.H"
#include "BvpOde.H"

double model_prob1_rhs( double x ){ return -1.0; }

int main( int argc, char* argv[] )
{

    SecondOrderODE
    ode_mp1
    (
        1.0,
        0.0,
        0.0,
        model_prob1_rhs,
        0.0,
        1.0
    );

    BoundaryConditions bc_mp1;
    bc_mp1.setLhsDirichletBc(0.0);
    bc_mp1.setRhsDirichletBc(0.0);

    BvpOde bvpode_mp1( &ode_mp1, &bc_mp1, 101 );
    bvpode_mp1.Solve();

    return 0;
}
