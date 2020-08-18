#include <iostream>
#include <cmath>
#include <string>

#include "IOTools.H"
#include "BvpODE.H"

double model_prob1_rhs(double x){return 1.0;}
double model_prob2_rhs(double x){return 34.0*sin(x);}

int main ( int argc, char* argv[] ){

    ScndOrderODE ode_mp1(-1.0, 0.0, 0.0, model_prob1_rhs, 0.0, 1.0);
    BoundaryConditions bc_mp1;
    bc_mp1.SetLHSDirichletBC(0.0);
    bc_mp1.SetRHSDirichletBC(0.0);

    BvpODE bvpode_mp1(&ode_mp1, &bc_mp1, 501);
    bvpode_mp1.SetFileName("model_problem_1.dat");
    bvpode_mp1.Solve();

    TecWriter tec;

    /*vars.push_back("X");
    vars.push_back("Y");
    vars.push_back("Z");

    tec.SetNumberOfVariables(5);
    tec.SetTitle("Testing object orientation");
    tec.SetVariables(vars);
    tec.WriteSolution<double>(d);
    tec.MakeHeader();*/

    return 0;
}
