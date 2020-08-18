#ifndef LINSYS_H
#define LINSYS_H

#include <cmath>
#include <iostream>
#include <fstream>
#include "Matrix.H"
#include "Vector.H"

class LinSys
{
private:
    Matrix* mpA;
    Vector* mpb;
    int mSize;
public:
    LinSys(Matrix& A, Vector& b);
    ~LinSys();

    virtual Vector Solve();
};

#endif
