#include <cassert>
#include "Matrix.H"

Matrix::Matrix(const Matrix& otherMatrix)
{
    mNumRows = otherMatrix.GetNumRows();
    mNumCols = otherMatrix.GetNumCols();
    mData = new double* [mNumRows];

    // allocating memory for mxn matrix
    for (int i = 0; i < mNumCols; i++)
    {
        mData[i] = new double [mNumCols];
    }

    for (int i = 0; i < mNumRows; i++)
    {
        for (int j = 0; j < mNumCols; j++)
        {
            mData[i][j] = otherMatrix.mData[i][j];
        }
    }
}

Matrix::Matrix(int numRows, int numCols)
{
    assert( numRows > 0 && numCols > 0 );
    mNumRows = numRows;
    mNumCols = numCols;
    mData = new double* [mNumRows];
    for (int i= 0; i < mNumRows; i++)
    {
        mData[i] = new double [mNumCols];
    }

    for (int i = 0; i < mNumRows; i++)
    {
        for (int j = 0; j < mNumCols; j++)
        {
            mData[i][j] = 0.0;
        }
    }
}

Matrix::~Matrix()
{
    // free memory allocated for object instantiation
    //std::cout << "Deleting Matrix object\n";

    for (int i = 0; i < mNumRows; i++)
    {
        delete[] mData[i];
    }

    delete[] mData;
}

int Matrix::GetNumRows() const
{
    return mNumRows;
}

int Matrix::GetNumCols() const
{
    return mNumCols;
}

double Matrix::Read(int i, int j) const
{
    return mData[i][j];
}

Matrix Matrix::Transpose()
{
    /* returns the trasnpose of a matrix to the current
    instance of the object. Given a matrix A with entries
    a_ij, that matrix will be turn into a matrix B such that
     b_ji = a_ij*/

    Matrix M(mNumCols,mNumRows);

    // attribute values for matrix elements
    for (int i = 0; i < mNumRows; i++)
    {
        for (int j = 0; j < mNumCols; j++)
        {
            M(j,i) = mData[i][j];

        }
    }
    return M;
}


double& Matrix::operator()(int i, int j)
{
    assert( i > -1 );
    assert( j > -1 );
    assert( i < mNumRows );
    assert( j < mNumCols );
    return mData[i][j];
}

bool Matrix::operator== (const Matrix& otherMatrix)
{
    assert( mNumRows == otherMatrix.mNumRows && mNumCols == otherMatrix.mNumCols );
    bool flag = true;

    for (int i = 0; i < mNumRows; i++)
    {
        for (int j = 0; j < mNumCols; j++)
        {
            if (mData[i][j] != otherMatrix.mData[i][j])
            {
                flag = false;
            }
        }
    }

    return flag;
}

Matrix& Matrix::operator= (const Matrix& otherMatrix)
{
    assert( mNumRows == otherMatrix.mNumRows && mNumCols == otherMatrix.mNumCols );

    for (int i = 0; i < mNumRows; i++)
    {
        for (int j = 0; j < mNumCols; j++)
        {
            mData[i][j] = otherMatrix.mData[i][j];
        }
    }

    return *this;
}

Matrix Matrix::operator+ () const
{
    Matrix M(mNumRows,mNumCols);

    for (int i = 0; i < mNumRows; i++)
    {
        for (int j = 0; j < mNumCols; j++)
        {
            M(i,j) = mData[i][j];
        }
    }

    return M;
}

Matrix Matrix::operator- () const
{
    Matrix M(mNumRows,mNumCols);

    for (int i = 0; i < mNumRows; i++)
    {
        for (int j = 0; j < mNumCols; j++)
        {
            M(i,j) = -mData[i][j];
        }
    }

    return M;
}

Matrix Matrix::operator+ (const Matrix& otherMatrix) const
{
    assert( mNumRows == otherMatrix.mNumRows && mNumCols == otherMatrix.mNumCols );
    Matrix M(mNumRows,mNumCols);

    for (int i = 0; i < mNumRows; i++)
    {
        for (int j = 0; j < mNumCols; j++)
        {
            M(i,j) = mData[i][j] + otherMatrix.mData[i][j];
        }
    }

    return M;
}

Matrix Matrix::operator- (const Matrix& otherMatrix) const
{
    assert( mNumRows == otherMatrix.mNumRows && mNumCols == otherMatrix.mNumCols );
    Matrix M(mNumRows,mNumCols);

    for (int i = 0; i < mNumRows; i++)
    {
        for (int j = 0; j < mNumCols; j++)
        {
            M(i,j) = mData[i][j] - otherMatrix.mData[i][j];
        }
    }

    return M;
}

Matrix Matrix::operator* (double a)
{
    Matrix M(mNumRows,mNumCols);

    // overloading * operator for scalar product
    for (int i = 0; mNumRows; i++)
    {
        for (int j = 0; j < mNumCols; j++)
        {
            M(i,j) = a * mData[i][j];
        }
    }
    return M;
}

Matrix Matrix::operator* (const Matrix& otherMatrix) const
{
    assert( mNumCols == otherMatrix.mNumRows );
    Matrix M(mNumRows,otherMatrix.mNumCols);
    double sum;

    for (int i = 0; i < mNumRows; i++)
    {
        for (int j = 0; j < otherMatrix.mNumCols; j++)
        {
            sum = 0.0;
            for (int k = 0; k < mNumCols; k++)
            {
                //M(i,j) += mData[i][k] * otherMatrix.mData[k][j];
                sum += mData[i][k] * otherMatrix.mData[k][j];
            }
            M(i,j) = sum;
        }
    }

    return M;
}

double Matrix::CalcDeterminant() const
{
    /* Computes the determinant of any NxN matrix using a
       recursive algorithm.*/
    assert( mNumRows == mNumCols );

    double determinant = 0.0;

    if (mNumRows == 1)
    {
        determinant = mData[0][0];
    }

    else
    {
        for (int i_outer = 0;  i_outer < mNumRows; i_outer++)
        {
            // creates a submatrix removing line and column
            Matrix m1(mNumRows-1, mNumCols-1);
            for (int i = 0; i < mNumRows-1; i++)
            {
                for (int j = 0; j < i_outer; j++)
                {
                    m1(i,j) = mData[i+1][j];
                }
                for (int k = i_outer; k < mNumCols-1; k++)
                {
                    m1(i,k) = mData[i+1][k+1];
                }
            }
            // recursively calls the determinant of the
            // submatrix
            double m1_determinant = m1.CalcDeterminant();
            determinant += pow(-1,i_outer) * mData[0][i_outer] * m1_determinant;
        }
    }

    return determinant;

}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Vector operator* (const Matrix& M, const Vector& v)
{
    assert( M.GetNumCols() == v.GetSize() );
    Vector x(M.GetNumRows());

    for (int i = 0; i < M.GetNumRows(); i++)
    {
        for (int j = 0; j < M.GetNumCols(); j++)
        {
            x(i) += M.Read(i,j) * v.Read(j);
        }
    }
    return x;
}

Vector operator* (const Vector& v, const Matrix& M)
{
    assert( M.GetNumRows() == v.GetSize() );
    Vector x(M.GetNumCols());

    for (int i = 0; i < M.GetNumCols(); i++)
    {
        for (int j = 0; j < M.GetNumRows(); j++)
        {
            x(i) += M.Read(j,i) * v.Read(j);
        }
    }
    return x;
}

void print(const Matrix& M)
{
    std::cout << "[ ";

    for (int i = 0; i < M.mNumRows; i++)
    {
        std::cout << "[ ";

        for (int j = 0; j < M.mNumCols; j++)
        {
            std::cout << M.mData[i][j] << " ";
        }

        std::cout << "] ";
    }

    std::cout << "]\n";
}

