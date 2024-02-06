#include <iostream>
#include <chrono>
#include <string>

#include "IOobject.h"
#include "BvpOde.h"

double model_prob1_rhs( double x ){ return 1.0; }

int main( int argc, char* argv[] )
{
  std::string SETUP_FILE = "../setup.inp";
  int imax = readInput<int>( "IMAX", SETUP_FILE );
  std::cout << "IMAX: " << imax << "\n";
  
  SecondOrderODE
  ode_mp1
  (
      0.0,                // coeff of Uxx
      1.0,                // coeff of Ux
      0.0,                // coeff of U
      model_prob1_rhs,    // pointer function to second order ODE
      0.0,                // xMin
      1.0                  // xMax
  );
  
  BoundaryConditions bc_mp1;
  bc_mp1.setLhsDirichletBc(1.0);
  bc_mp1.setRhsDirichletBc(0.0);
  
  BvpOde bvpode_mp1( &ode_mp1, &bc_mp1, imax );
  bvpode_mp1.Solve();
  
  return 0;
}
