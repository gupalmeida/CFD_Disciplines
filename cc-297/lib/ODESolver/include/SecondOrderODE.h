#ifndef SECONDORDERODE_H
#define SECONDORDERODE_H

#include <iostream>
#include <cmath>
#include <vector>

#include "Vector.h"
#include "Matrix.h"

class SecondOrderODE
{
    // Boundary value class is able to access
    // the coefficients, etc of this equation
    friend class BvpOde;

    private:
        // coefficients on LHS of ODE
        double mCoeffOfUxx;
        double mCoeffOfUx;
        double mCoeffOfU;

        // function on RHS of ODE
        double (*mpRhsFunction)(double x);

        // interval of the domain
        double mXmin;
        double mXmax;

    public:
        // Constructor
        SecondOrderODE
        (
            double coeffUxx,
            double coeffUx,
            double coeffU,
            double (*rightHandSide)(double),
            double xMin,
            double xMax
        )
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
