#include <cassert>
#include "BoundaryConditions.H"

BoundaryConditions::BoundaryConditions()
{
    mLHSDirichlet = false;
    mRHSDirichlet = false;
    mLHSNeumann = false;
    mRHSNeumann = false;
}

void BoundaryConditions::SetLHSDirichletBC(double lhsValue)
{
    assert(!mLHSDirichlet);
    mLHSDirichlet = true;
    mLHSBCValue = lhsValue;
}

void BoundaryConditions::SetRHSDirichletBC(double rhsValue)
{
    assert(!mRHSDirichlet);
    mRHSDirichlet = true;
    mRHSBCValue = rhsValue;
}

void BoundaryConditions::SetLHSNeumannBC(double lhsDiffValue)
{
    assert(!mLHSNeumann);
    mLHSNeumann = true;
    mLHSBCValue = lhsDiffValue;
}

void BoundaryConditions::SetRHSNeumannBC(double rhsDiffValue)
{
    assert(!mRHSNeumann);
    mRHSNeumann = true;
    mRHSBCValue = rhsDiffValue;
}

