#ifndef BVPODE_H
#define BVPODE_H

#include <cmath>
#include <string>

#include "IOTools.H"
#include "Matrix.H"
#include "Vector.H"
#include "LinSys.H"
#include "FDGrid.H"
#include "BoundaryConditions.H"
#include "ScndOrderODE.H"

class BvpODE
{
private:
    BvpODE(const BvpODE& otherBvpODE);

    int mNumNodes;
    std::string mFileName;

    FDGrid* mpGrid;
    ScndOrderODE* mpODE;
    BoundaryConditions* mpBCs;
    Vector* mpSolution;
    Vector* mpRHSVec;
    Matrix* mpLHSMat;
    LinSys* mpLinSys;

    void PopulateVector();
    void PopulateMatrix();
    void ApplyBCs();

public:
    BvpODE(ScndOrderODE* pODE, BoundaryConditions* pBCs,
                                        int numNodes);
    ~BvpODE();
    
    void SetFileName(const std::string name)
    {
        mFileName = name;
    }

    void Solve();
    //void WriteSolutionFile();
};

#endif
