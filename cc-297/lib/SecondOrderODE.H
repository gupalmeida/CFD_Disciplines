#ifndef SECONDORDERODE_H
#define SECONDORDERODE_H

#include <iostream>
#include <cmath>
#include <vector>

#include "Vector.H"
#include "Matrix.H"

class SecondOrderODE
{
    private:
        // coefficients on LHS of ODE
        double mCoeffOfUxx;
        double mCoeffOfUx;
        double mCoeffOfU;

        // function on RHS of ODE
        double (*mpRhsFunction)(double x);

        // interval of the domain
        double mXmin;
        double mXMax;

    public:
        SecondOrderODE( double coeffUxx, double coeffUx,
                        double coeffU,
                        double (*rightHandSide)(double),
                        double xMin, double xMax)
        {
            mXmin = xMin;
            mXmax = xMax;
            mpRhsFunction = rightHandSide;
            mCoeffOfUxx = coeffUxx;
            mCoeffOfUx = coeffUx;
            mCoeffOfU = coeffU;
        }
};

#endif
