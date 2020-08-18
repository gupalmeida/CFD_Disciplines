#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <vector>
#include <cmath>
#include "Vector.H"

class Matrix
{
private:
    double** mData;
    int mNumRows;
    int mNumCols;

public:
    Matrix(const Matrix& otherMatrix);
    Matrix(int numRows, int numCols);
    ~Matrix();
    int GetNumRows() const;
    int GetNumCols() const;
    double Read(int i, int j) const;
    double CalcDeterminant() const;
    Matrix Transpose();
    double& operator()(int i,int j);
    bool operator== (const Matrix& otherMatrix);
    Matrix& operator= (const Matrix& otherMatrix);
    
    // unary operations
    Matrix operator+ () const;
    Matrix operator- () const;

    // binary operations
    Matrix operator+ (const Matrix& otherMatrix) const;
    Matrix operator- (const Matrix& otherMatrix) const;
    Matrix operator* (double a);
    Matrix operator* (const Matrix& otherMatrix) const;

    friend Vector operator* (const Matrix& M, const Vector& v);
    friend Vector operator* (const Vector& v, const Matrix& M);
    friend void print(const Matrix& M);
};

#endif
