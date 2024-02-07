#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <cmath>
#include <cassert>
using namespace std;

#include "Vector.h"

class Matrix
{
private:
  double** mData; // matrix data
  int mRows, mCols; // matrix dimensions MxN

public:
  Matrix( int numRows, int numCols ); // new matrix constructor
  Matrix( const Matrix& otherMatrix ); // new matrix constructor
  ~Matrix(); // matrix destructor
  // getting and attributing data
  int getNumRows() const;
  int getNumCols() const;
  void ones();
  void eye();
  void setUpperDiag( double u=1.0 );
  void setDiag( double d=1.0 );
  void setLowerDiag( double l=1.0 );
  // assignment
  //double& operator[](int i, int j); // zero-based indexing
  double& operator()(int i, int j); // zero-based indexing
  Matrix& operator=( const Matrix& otherMatrix );
  Matrix operator+(  );
  Matrix operator-(  );
  Matrix operator+( const Matrix& otherMatrix );
  Matrix operator-( const Matrix& otherMatrix );
  Matrix operator*( double a );

  // matrix operations
  Matrix T(  ); // returns the transposed matrix
  //Matrix& dot( const Matrix& n );
  double Det(  ) const;

  friend void print( const Matrix& m );
};

// python-like friend function to print out matrix values
void print( const Matrix& m );
void print( const double matrixEntry );

#endif
