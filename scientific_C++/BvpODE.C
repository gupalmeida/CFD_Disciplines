#include <cassert>
#include "BvpODE.H"

BvpODE::BvpODE(ScndOrderODE* pODE, BoundaryConditions* pBCs,
                                        int numNodes)
{
    mpODE = pODE;
    mpBCs = pBCs;

    mNumNodes = numNodes;

    mpGrid = new FDGrid(mNumNodes, pODE->mXMin, pODE->mXMax);
    mpSolution = new Vector(mNumNodes);
    mpRHSVec = new Vector(mNumNodes);
    mpLHSMat = new Matrix(mNumNodes, mNumNodes);

    mFileName = "ode_output.dat";
    mpLinSys = NULL;
}

BvpODE::~BvpODE()
{
    delete mpSolution;
    delete mpRHSVec;
    delete mpLHSMat;
    delete mpGrid;
    if (mpLinSys)
    {
        delete mpLinSys;
    }
}

void BvpODE::Solve()
{
    PopulateVector();
    PopulateMatrix();
    ApplyBCs();
    mpLinSys = new LinSys(*mpLHSMat, *mpRHSVec);
    *mpSolution = mpLinSys->Solve();
    //WriteSolutionFile();
    print(*mpSolution);
}

void BvpODE::PopulateMatrix()
{
    for (int i = 1; i < mNumNodes-1; i++)
    {
        double xm = mpGrid->mNodes[i-1].xCoord;
        double x = mpGrid->mNodes[i].xCoord;
        double xp = mpGrid->mNodes[i+1].xCoord;
        double alpha = 2.0/((xp-xm)*(x-xm));
        double beta = -2.0/((xp-x)*(x-xm));
        double gamma = 2.0/((xp-xm)*(xp-x));

        (*mpLHSMat)(i,i-1) = (mpODE->mCoeffUxx)*alpha -
        (mpODE->mCoeffUx)/(xp-xm);
        (*mpLHSMat)(i,i) = (mpODE->mCoeffUxx)*beta +
            mpODE->mCoeffU;
        (*mpLHSMat)(i,i+1) = (mpODE->mCoeffUxx)*gamma +
            (mpODE->mCoeffUx)/(xp-xm);
    }
}

void BvpODE::PopulateVector()
{
    for (int i=1; i<mNumNodes-1; i++)
    {
        double x = mpGrid->mNodes[i].xCoord;
        (*mpRHSVec)(i) = mpODE->mpRHS(x);
    }
}

void BvpODE::ApplyBCs()
{
    bool left_bc_applied = false;
    bool right_bc_applied = false;

    if (mpBCs->mLHSDirichlet)
    {
        (*mpLHSMat)(0,0) = 1.0;
        (*mpRHSVec)(0) = mpBCs->mLHSBCValue;
        left_bc_applied = true;
    }

    if (mpBCs->mRHSDirichlet)
    {
        (*mpLHSMat)(mNumNodes-1,mNumNodes-1) = 1.0;
        (*mpRHSVec)(mNumNodes-1) = mpBCs->mRHSBCValue;
        right_bc_applied = true;
    }

    if (mpBCs->mLHSNeumann)
    {
        assert( left_bc_applied == false );

        double h = mpGrid->mNodes[1].xCoord - mpGrid->mNodes[0].xCoord;

        (*mpLHSMat)(0,0) = -1.0/h;
        (*mpLHSMat)(0,1) = 1.0/h;
        (*mpRHSVec)(0) = mpBCs->mLHSBCValue;
        left_bc_applied = true;
    }

    if (mpBCs->mRHSNeumann)
    {
        assert( right_bc_applied == false );

        double h = mpGrid->mNodes[mNumNodes-1].xCoord
                            - mpGrid->mNodes[mNumNodes-2].xCoord;

        (*mpLHSMat)(mNumNodes-1,mNumNodes-2) = -1.0/h;
        (*mpLHSMat)(mNumNodes-1,mNumNodes-1) = 1.0/h;
        (*mpRHSVec)(mNumNodes-1) = mpBCs->mRHSBCValue;
        right_bc_applied = true;
    }

    assert( left_bc_applied );
    assert( right_bc_applied );
}
