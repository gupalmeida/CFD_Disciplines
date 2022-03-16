#include <iostream>
#include <cmath>
#include <cassert>
#include "Matrix.H"

// new matrix constructor
Matrix::Matrix( int numRows, int numCols )
{
    assert( numRows > 0 );
    assert( numCols > 0 );
    mRows = numRows;
    mCols = numCols;
    mData = new double* [mRows];
    for ( int i = 0; i <  mRows; i++)
    {
        mData[i] = new double [mCols]; // allocates memory for new matrix
    }
    for (int i = 0; i < mRows; i++)
    {
        for (int j = 0; j < mCols; j++)
        {
            mData[i][j] = 0.0; // initializes new matrix with zeros
        }
    }
}

Matrix::Matrix( const Matrix& otherMatrix )
{
    assert( otherMatrix.mRows > 0 );
    assert( otherMatrix.mCols > 0 );
    mRows = otherMatrix.mRows;
    mCols = otherMatrix.mCols;
    mData = new double* [mRows];
    for ( int i = 0; i <  mRows; i++)
    {
        mData[i] = new double [mCols]; // allocates memory for new matrix
    }
    for (int i = 0; i < mRows; i++)
    {
        for (int j = 0; j < mCols; j++)
        {
            mData[i][j] = otherMatrix.mData[i][j]; // copy entries from otherMatrix
        }
    }
}

// destructor
Matrix::~Matrix()
{
    for ( int i = 0; i < mRows; i++ )
    {
        delete[] mData[i];
    }
    delete[] mData;
}

// getting and attributing data
int Matrix::getNumRows(  ) const
{
    return mRows;
}

int Matrix::getNumCols(  ) const
{
    return mCols;
}

void Matrix::ones(  )
{
    // set all matrix entries to 1.0
    assert( mRows > 0 );
    assert( mCols > 0 );
    for ( int i = 0; i < mRows; i++ )
    {
        for ( int j = 0; j < mCols; j++ )
        {
            mData[i][j] = 1.0;
        }
    }
}

void Matrix::eye(  )
{
    // set to identity matrix
    assert( mRows > 0 );
    assert( mCols > 0 );
    // first turn all entries to 0
    for ( int i = 0; i < mRows; i++ )
    {
        for ( int j = 0; j < mCols; j++ )
        {
            mData[i][j] = 0.0;
        }
    }
    // then changes diagonal values to 1.0
    for ( int i = 0; i < mRows; i++ )
    {
        for ( int j = i; j < i+1; j++ )
        {
            mData[i][j] = 1.0;
        }
    }
}

void Matrix::setUpperDiag( double u )
{
    // set upper diagonal to user defined value.
    // default is 1.0.
    assert( mRows > 0 );
    assert( mCols > 0 );
    assert( mRows == mCols);
    for ( int i = 0; i < mRows-1; i++ )
    {
        for ( int j = i; j < i+1; j++ )
        {
            mData[i][j+1] = u;
        }
    }
}

void Matrix::setDiag( double d )
{
    // set diagonal to user defined value.
    // default is 1.0.
    assert( mRows > 0 );
    assert( mCols > 0 );
    assert( mRows == mCols );
    for ( int i = 0; i < mRows; i++ )
    {
        for ( int j = i; j < i+1; j++ )
        {
            mData[i][j] = d;
        }
    }
}

void Matrix::setLowerDiag( double l )
{
    // set lower diagonal to user defined value.
    // default is 1.0.
    assert( mRows > 0 );
    assert( mCols > 0 );
    assert( mRows == mCols );
    for ( int i = 1; i < mRows; i++ )
    {
        for ( int j = i; j < i+1; j++ )
        {
            mData[i][j-1] = l;
        }
    }
}


// assignment
//double& Matrix::operator[]( int i, int j )
//{
// operator [] does not accept two entries like i,j.
//    assert( i > -1 );
//    assert( i < mRows );
//    assert( j > -1 );
//    assert( j < mCols );
//    return mData[i][j];
//}

double& Matrix::operator()( int i, int j )
{
    assert( i > -1 );
    assert( i < mRows );
    assert( j > -1 );
    assert( j < mCols );
    return mData[i][j];
}

Matrix& Matrix::operator=( const Matrix& otherMatrix )
{
    assert( mRows == otherMatrix.mRows );
    assert( mCols == otherMatrix.mCols );
    for ( int i = 0; i < mRows; i++ )
    {
        for ( int j = 0; j < mCols; j++ )
        {
            mData[i][j] = otherMatrix.mData[i][j];
        }
    }
    return *this;
}

Matrix Matrix::operator+( )
{
    assert( mRows > 0 );
    assert( mCols > 0 );
    Matrix m(mRows,mCols);
    for ( int i = 0; i < mRows; i++ )
    {
        for ( int j = 0; j < mCols; j++ )
        {
            m(i,j) = mData[i][j];
        }
    }
    return m;
}

Matrix Matrix::operator-( )
{
    assert( mRows > 0 );
    assert( mCols > 0 );
    Matrix m(mRows,mCols);
    for ( int i = 0; i < mRows; i++ )
    {
        for ( int j = 0; j < mCols; j++ )
        {
            m(i,j) = - mData[i][j];
        }
    }
    return m;
}

Matrix Matrix::operator+( const Matrix& otherMatrix )
{
    assert( mRows == otherMatrix.mRows );
    assert( mCols == otherMatrix.mCols );
    Matrix m(mRows,mCols);
    for ( int i = 0; i < mRows; i++ )
    {
        for ( int j = 0; j < mCols; j++ )
        {
            m(i,j) = mData[i][j] + otherMatrix.mData[i][j];
        }
    }
    return m;
}

Matrix Matrix::operator-( const Matrix& otherMatrix )
{
    assert( mRows == otherMatrix.mRows );
    assert( mCols == otherMatrix.mCols );
    Matrix m(mRows,mCols);
    for ( int i = 0; i < mRows; i++ )
    {
        for ( int j = 0; j < mCols; j++ )
        {
            m(i,j) = mData[i][j] - otherMatrix.mData[i][j];
        }
    }
    return m;
}

Matrix Matrix::operator*( double a )
{
    assert( mRows == mCols );
    Matrix m(mRows,mCols);
    for ( int i = 0; i < mRows; i++ )
    {
        for ( int j = 0; j < mCols; j++ )
        {
            m(i,j) = a * mData[i][j];
        }
    }
    return m;
}

// matrix operations
Matrix Matrix::T(  )
{
    assert( mRows > 0 );
    assert( mCols > 0 );
    Matrix m(mCols,mRows);
    for ( int i = 0; i < mRows; i++ )
    {
        for ( int j = 0; j < mCols; j++ )
        {
            m(j,i) = mData[i][j];
        }
    }
    return m;
}

double Matrix::Det(  ) const
{
    assert( mRows == mCols );
    double determinant = 0.0;
    if ( mRows == 1 )
    {
        determinant = mData[0][0];
        return determinant;
    }
    else
    {
        for (int i_outer = 0;  i_outer < mRows; i_outer++)
        {
            Matrix m1(mRows-1, mCols-1);
            for (int i = 0; i < mRows-1; i++)
            {
                for (int j = 0; j < i_outer; j++)
                {
                    m1(i,j) = mData[i+1][j];
                }
                for (int k = i_outer; k < mCols-1; k++)
                {
                    m1(i,k) = mData[i+1][k+1];
                }
            }
                double m1_determinant = m1.Det();
                determinant += pow(-1,i_outer) * mData[0][i_outer] * m1_determinant;
        }
    }
        return determinant;
}

//Matrix& Matrix::dot( const Matrix& n )
//{
//    assert( mCols == n.mRows );
//    Matrix r(mRows,n.mCols);
//    double sum;
//
//    for (int i = 0; i < mRows; i++)
//    {
//        for (int j = 0; j < n.mCols; j++)
//        {
//            sum = 0.0;
//
//            for (int k = 0; k < mCols; k++)
//            {
//                sum += (*this)(i,k) * n(k,j);
//            }
//
//            r(i,j) = sum;
//        }
//    }
//
//    (*this) = r;
//
//    return *this;
//}
// python-like friend function to print out matrix values
void print( const Matrix& m )
{
    assert( m.mRows > 0 );
    assert( m.mCols > 0 );
    for ( int i = 0; i < m.mRows; i++ )
    {
        for ( int j = 0; j < m.mCols; j++ )
        {
            std::cout << m.mData[i][j] << '\t';
        }
        std::cout << '\n';
    }
}

void print( const double matrixEntry )
{
    // prints out a single entry of a
    // matrix
    std::cout << matrixEntry << '\n';
}

