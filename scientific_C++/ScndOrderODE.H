#ifndef SCNDORDERODE_H
#define SCNDORDERODE_H

#include <cmath>

class ScndOrderODE
{
public:
    friend class BvpODE;

private:
    double mCoeffUxx;
    double mCoeffUx;
    double mCoeffU;

    double (*mpRHS) (double x);
    
    double mXMax;
    double mXMin;

public:
    ScndOrderODE(double coeffUxx, double coeffUx,
                    double coeffU,
                    double (*rightHandSide)(double),
                    double xMin, double xMax)
    {
        mCoeffUxx = coeffUxx;
        mCoeffUx = coeffUx;
        mCoeffU = coeffU;
        mpRHS = rightHandSide;
        mXMin = xMin;
        mXMax = xMax;
    }
};

#endif
