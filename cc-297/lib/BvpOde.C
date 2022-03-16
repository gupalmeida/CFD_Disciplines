#include "BvpOde.H"


BvpOde::BvpOde(SecondOrderOde* pOde, BoundaryConditions* pBcs, int numNodes)
{

    mpOde = pOde;

    mpBCs = pBcs;

    mNumOfNodes = numNodes;

    mpMesh = new Mesh( pOde->xMin, pOde->xMax, mNumOfNodes );
    mpSolution = new Vector( mNumOfNodes );
    mpRHS = new Vector( mNumOfNodes );
    mpLinMat = new Matrix( mNumOfNodes, mNumOfNodes );
    mpLinSys = NULL;
}

BvpOde::~BvpOde()
{
    
}
