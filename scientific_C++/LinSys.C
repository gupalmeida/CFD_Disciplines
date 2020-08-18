#include <cassert>
#include "LinSys.H"

LinSys::LinSys(Matrix& A, Vector& b)
{
    assert( A.GetNumCols() == A.GetNumRows() );
    assert( A.GetNumRows() == b.GetSize() );
    
    mSize = A.GetNumRows();
    // creates local copies to A and b to avoid changing
    // their original values during the process of solving
    // system of linear equations
    mpA = new Matrix(A);
    mpb = new Vector(b);

}

LinSys::~LinSys()
{
    delete mpA;
    delete mpb;
}

Vector LinSys::Solve()
{
    /* Solves the linear system of equations using
       Gauss Elimination without pivoting */

    Vector m( mSize ); // vector for storing forward sweep coeffs
    Vector x( mSize ); // solution vector

    // referencing matrix and vector to make code readable
    Matrix& rA = *mpA;
    Vector& rb = *mpb;

    // forward sweep
    for (int k = 0; k < mSize-1; k++)
    {
        for (int i = k+1; i < mSize; i++)
        {
            m(i) = rA(i,k) / rA(k,k);
            for (int j = k; j < mSize; j++)
            {
                rA(i,j) -= rA(k,j) * m(i);
            }
            rb(i) -= rb(k) * m(i);
        }
    }

    // backward substitution
    for (int i = mSize-1; i > -1; i--)
    {
        x(i) = rb(i);
        for (int j = i+1; j < mSize; j++)
        {
            x(i) -= rA(i,j) * x(j);
        }
        x(i) /= rA(i,i);
    }

    return x;
}
