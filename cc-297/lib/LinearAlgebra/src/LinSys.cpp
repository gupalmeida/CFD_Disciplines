#include "LinSys.h"


LinSys::LinSys(const Matrix& A, const Vector& b)
{
  // check matrix and vector are of same size
  assert( A.getNumCols() == mSize );
  assert( b.size() == mSize );

  // define size of solution vector same
  // as associated matrix
  mSize = A.getNumRows();

  mpA = new Matrix(A); // using copy constructor
  mpb = new Vector(b); // using copy constructor
}

LinSys::~LinSys()
{
  delete mpA;
  delete mpb;
}

Vector LinSys::Solve()
{
  // The solution method adopted in the Linear System class
  // is the Gauss elimination without pivoting as default.
  // Other linear solvers are further implemented in other
  // classes.

  // creates and initializes solution vector
  Vector x(mSize);

  // creates coefficient vector for Gauss elimination procedure
  Vector m(mSize);

  // getting references for the passed Matrix and Vector
  Matrix rA = *mpA;
  Vector rb = *mpb;

  // forward sweep
  for (int k = 0 ; k < mSize-1; k++ )
  {
      for ( int i = k+1; i < mSize; i++ )
      {
          m(i) = rA(i,k)/rA(k,k);
          for ( int j=k ; j < mSize; j++)
          {
              rA(i,j) -= m(i) * rA(k,j);
          }
          rb(i) -= m(i) * rb(k);
      }
  }

  // backward substitution
  x(mSize-1) = rb(mSize-1)/rA(mSize-1,mSize-1);
  for ( int i = mSize-2; i > -1; i-- )
  {
      x(i) = rb(i);
      for ( int j=i+1; j < mSize; j++ )
      {
          x(i) -= rA(i,j) * x(j);
      }
      x(i) /= rA(i,i);
  }

  return x;
}
