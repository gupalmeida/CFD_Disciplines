#ifndef BOUNDARYCONDITIONS_H
#define BOUNDARYCONDITIONS_H

#include <cmath>

class BoundaryConditions
{
public:
    friend class BvpODE;

private:
    bool mLHSDirichlet;
    bool mRHSDirichlet;
    bool mLHSNeumann;
    bool mRHSNeumann;
    double mLHSBCValue;
    double mRHSBCValue;

public:
    BoundaryConditions();
    void SetLHSDirichletBC(double lhsValue);
    void SetRHSDirichletBC(double rhsValue);
    void SetLHSNeumannBC(double lhsDiffValue);
    void SetRHSNeumannBC(double rhsDiffValue);
};

#endif
