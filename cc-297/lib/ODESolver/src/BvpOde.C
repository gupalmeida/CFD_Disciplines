#include "BvpOde.H"


BvpOde::BvpOde
(
    SecondOrderODE* pOde,
    BoundaryConditions* pBcs,
    int numNodes
)
{
    mpOde = pOde;
    mpBCs = pBcs;

    mNumOfNodes = numNodes;

    mpMesh = new Mesh( pOde->mXmin, pOde->mXmax, mNumOfNodes );
    mpSolution = new Vector( mNumOfNodes );
    mpRHS = new Vector( mNumOfNodes );
    mpLinMat = new Matrix( mNumOfNodes, mNumOfNodes );
    mpLinSys = NULL;
}

BvpOde::~BvpOde()
{
    // free allocated memory
    delete mpSolution;
    delete mpRHS;
    delete mpLinMat;
    delete mpMesh;
    // deletes only if Solve has been called
    if (mpLinSys)
    {
        delete mpLinSys;
    }
}

void BvpOde::Solve()
{
    populateMatrix();
    populateVector();
    applyBCs();
    
    mpLinSys = new LinSys( *mpLinMat, *mpRHS );
    *mpSolution = mpLinSys->Solve();
    print( *mpSolution );
    writeSolutionToFile( "solution.out" );
}

void BvpOde::populateMatrix()
{
    // loop trhough all internal points of the
    // computational mesh
    for ( int i = 1; i < mNumOfNodes-1; i++ )
    {
        // xm, x and xp are, respectively, x-1, x and x+1
        double xm = mpMesh->mNodes[i-1].xCoord;
        double x = mpMesh->mNodes[i].xCoord;
        double xp = mpMesh->mNodes[i+1].xCoord;

        // coefficients using central differencing scheme
        // alpha * x_{i-1} + beta * x + gamma * x_{i+1}
        double alpha = 2.0 / ((xp-xm)*(x-xm));
        double beta = -2.0 / ((xp-x)*(x-xm));
        double gamma = 2.0 / ((xp-xm)*(xp-x));

        (*mpLinMat)(i,i-1) = (mpOde->mCoeffOfUxx)*alpha -
            (mpOde->mCoeffOfUx) / (xp-xm);
        (*mpLinMat)(i,i) = (mpOde->mCoeffOfUxx)*beta +
            (mpOde->mCoeffOfU);
        (*mpLinMat)(i,i+1) = (mpOde->mCoeffOfUxx)*gamma +
            (mpOde->mCoeffOfUx) / (xp-xm);
    }

}

void BvpOde::populateVector()
{
    // loop trhough elements corresponding to internal
    // points of computational mesh
    for (int i = 1; i < mNumOfNodes-1; i++)
    {
        double x = mpMesh->mNodes[i].xCoord;
        (*mpRHS)(i) = mpOde->mpRhsFunction(x);
    }
}

void BvpOde::applyBCs()
{
    bool left_bc_applied = false;
    bool right_bc_applied = false;

    if (mpBCs->mLhsBcIsDirichlet)
    {
        (*mpLinMat)(0,0) = 1.0;
        (*mpRHS)(0) = mpBCs->mLhsBcValue;
        left_bc_applied = true;
    }

    if (mpBCs->mRhsBcIsDirichlet)
    {
        (*mpLinMat)(mNumOfNodes-1,mNumOfNodes-1) = 1.0;
        (*mpRHS)(mNumOfNodes-1) = mpBCs->mRhsBcValue;
        right_bc_applied = true;
    }

    if (mpBCs->mLhsBcIsNeumann)
    {
        assert( left_bc_applied == false );

        double h = mpMesh->mNodes[1].xCoord - mpMesh->mNodes[0].xCoord;

        (*mpLinMat)(0,0) = -1.0/h;
        (*mpLinMat)(0,1) = 1.0/h;
        (*mpRHS)(0) = mpBCs->mLhsBcValue;
        left_bc_applied = true;
    }

    if (mpBCs->mRhsBcIsNeumann)
    {
        assert( right_bc_applied == false );

        double h = mpMesh->mNodes[mNumOfNodes-1].xCoord - mpMesh->mNodes[mNumOfNodes-2].xCoord;

        (*mpLinMat)(mNumOfNodes-1,mNumOfNodes-2) = -1.0/h;
        (*mpLinMat)(mNumOfNodes-1,mNumOfNodes-1) = 1.0/h;
        (*mpRHS)(mNumOfNodes-1) = mpBCs->mLhsBcValue;
        right_bc_applied = true;
    }

    assert( left_bc_applied = true );
    assert( right_bc_applied = true );
}

// friend functions
void BvpOde::writeSolutionToFile( string fname )
{
    ofstream outputFile( fname );

    for ( int i=0; i < mNumOfNodes; i++ )
    {
        outputFile<< (*mpSolution)(i) << "\n";
    }

    outputFile.close();
}
