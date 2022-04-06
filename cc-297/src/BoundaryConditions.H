#ifndef BOUNDARYCONDITIONS_H
#define BOUNDARYCONDITIONS_H

#include <iostream>
#include <cmath>
#include <cassert>

#include "Aux.H"

class BoundaryConditions
{
    public:
        friend class BvpOde;

    private:
        bool mLhsBcIsDirichlet;
        bool mLhsBcIsNeumann;
        bool mRhsBcIsDirichlet;
        bool mRhsBcIsNeumann;
        double mLhsBcValue;
        double mRhsBcValue;

    public:
        BoundaryConditions();

        void setLhsDirichletBc( double lhsBcValue );
        void setLhsNeumannBc( double lhsBcValue );
        void setRhsDirichletBc( double rhsBcValue );
        void setRhsNeumannBc( double rhsBcValue );
};

#endif
