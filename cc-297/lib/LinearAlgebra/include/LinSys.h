#ifndef LINSYS_H
#define LINSYS_H

#include <iostream>
#include <cmath>
#include <cassert>
using namespace std;

#include "Matrix.h"
#include "Vector.h"

class LinSys
{
    private:
        int mSize;
        Matrix* mpA;
        Vector* mpb;

    public:
        LinSys(const Matrix& A, const Vector& b);
        ~LinSys();

        virtual Vector Solve();
};

#endif
