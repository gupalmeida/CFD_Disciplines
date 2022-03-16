#ifndef LINSYS_H
#define LINSYS_H

#include "Matrix.H"
#include "Vector.H"

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
