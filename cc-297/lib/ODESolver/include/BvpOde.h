#ifndef BVPODE_H
#define BVPODE_H

#include <iostream>
#include <cmath>
#include <cassert>
#include <string>

#include "Vector.h"
#include "Matrix.h"
#include "LinSys.h"
#include "Mesh.h"
#include "SecondOrderODE.h"
#include "BoundaryConditions.h"
#include "IOobject.h"

class BvpOde
{
    private:
        // making copy constructor private for only allow construction
        // from PDE, boundary conditions and number of nodes in mesh.
        BvpOde( const BvpOde& otherBvpOde){}

        // number of nodes in the mesh and pointer to mesh
        int mNumOfNodes;
        Mesh* mpMesh;

        // pointer to an instance of an ODE
        SecondOrderODE* mpOde;

        // pointer to an instance of BC
        BoundaryConditions* mpBCs;

        // pointer to the solution vector
        Vector* mpSolution;

        // pointer to RHS
        Vector* mpRHS;

        // pointer to matrix of coefficients
        Matrix* mpLinMat;

        // pointer to the resulting linear system
        LinSys* mpLinSys;

        // methods for setting up linear system
        void populateMatrix();
        void populateVector();
        void applyBCs();

    public:
        BvpOde(SecondOrderODE* pOde, BoundaryConditions* pBcs, int numNodes);
        ~BvpOde();

        void Solve();
        void writeSolutionToFile( string fname );
};

#endif
