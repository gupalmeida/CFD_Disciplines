#include "BoundaryConditions.H"

BoundaryConditions::BoundaryConditions()
{
    mLhsBcIsDirichlet = false;
    mLhsBcIsNeumann   = false;
    mRhsBcIsDirichlet = false;
    mRhsBcIsNeumann   = false;

    mLhsBcValue = 0.0;
    mRhsBcValue = 0.0;
}

void BoundaryConditions::setLhsDirichletBc( double lhsBcValue )
{
    assert( !mLhsBcIsNeumann );
    mLhsBcIsDirichlet = true;
    mLhsBcValue = lhsBcValue;
}

void BoundaryConditions::setLhsNeumannBc( double lhsBcValue )
{
    assert( !mLhsBcIsDirichlet );
    mLhsBcIsNeumann = true;
    mLhsBcValue = lhsBcValue;
}

void BoundaryConditions::setRhsDirichletBc( double rhsBcValue )
{
    assert( !mRhsBcIsNeumann );
    mRhsBcIsDirichlet = true;
    mRhsBcValue = rhsBcValue;
}

void BoundaryConditions::setRhsNeumannBc( double rhsBcValue )
{
    assert( !mRhsBcIsDirichlet );
    mRhsBcIsNeumann = true;
    mRhsBcValue = rhsBcValue;
}

