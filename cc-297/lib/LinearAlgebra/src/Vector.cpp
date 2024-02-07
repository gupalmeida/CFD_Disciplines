#include "Vector.h"

// constructor for vector with given size
// initializing the new instance with all zeros
Vector::Vector( int size )
{
  assert( size > 0 );
  mSize = size;
  mData = new double [mSize];
  for (int i = 0; i < mSize; i++)
  {
    mData[i] = 0.0;
  }
}

// constructor for when copying another vector
Vector::Vector( const Vector& v )
{
  assert( v.mSize > 0 );
  mSize = v.mSize;
  mData = new double [mSize];
  for (int i = 0; i < mSize; i++)
  {
    mData[i] = v.mData[i];
  }
}

// destructor
Vector::~Vector()
{
  delete[] mData;
}

// vector initialization
void Vector::ones()
{
  assert( mSize > 0 );
  for ( int i = 0; i < mSize; i++ )
  {
    mData[i] = 1.0;
  }
}

// vector info and printing
int Vector::size() const
{
  return mSize;
}

void Vector::print() const
{
  assert( mSize > 0 );
  for (int i = 0; i < mSize; i++)
  {
    std::cout << mData[i] << '\n';
  }
}

void Vector::print( int i ) const
{
    // verify index is within bounds
    assert( i > -1 );
    assert( i < mSize );
    std::cout << mData[i] << '\n';
}

// reading vector element
double& Vector::operator[]( int i )
{
    // zero-based indexing
    // verify index is within bounds
    assert( i > -1 );
    assert( i < mSize );
    return mData[i];
}

double& Vector::operator()( int i )
{
    // zero-based indexing
    assert( i > -1 );
    assert( i < mSize );
    return mData[i];
}

// assignment and operations
Vector& Vector::operator=( const Vector& v )
{
    assert( mSize == v.mSize );
    for ( int i = 0; i < mSize; i++ )
    {
        mData[i] = v.mData[i];
    }
    return *this;
}

Vector Vector::operator+() const        // unary operator
{
    assert( mSize > 0 );
    Vector u(mSize);
    for ( int i = 0; i < mSize; i++ )
    {
        u[i] = mData[i];
    }
    return u;
}

Vector Vector::operator-() const        // unary operator
{
    assert( mSize > 0 );
    Vector u(mSize);
    for ( int i = 0; i < mSize; i++ )
    {
        u[i] = -mData[i];
    }
    return u;
}

Vector Vector::operator+( const Vector& v ) const // binary operator
{
    assert( mSize == v.mSize );
    Vector u(mSize);
    for ( int i = 0; i < mSize; i++ )
    {
        u[i] = mData[i] + v.mData[i];
    }
    return u;
}

Vector Vector::operator-( const Vector& v ) const // binary operator
{
    assert( mSize == v.mSize );
    Vector u(mSize);
    for ( int i = 0; i < mSize; i++ )
    {
        u[i] = mData[i] - v.mData[i];
    }
    return u;
}

Vector Vector::operator*( double a ) const
{
    Vector v(mSize);
    for ( int i = 0; i < mSize; i++ )
    {
        v[i] = a * mData[i];
    }
    return v;
}

double Vector::pNorm( int p ) const
{
    double norm, sum = 0.0;
    for ( int i = 0; i < mSize; i++ )
    {
        sum += pow( fabs( mData[i] ), p );
    }
    norm = pow( sum, (1.0/ ((double) p)) );
    return norm;
}

//Vector& Vector::operator*( double a )
//{
//    assert( mSize > 0 );
//    for ( int i = 0; i < mSize; i++ )
//    {
//        mData[i] = a * mData[i];
//    }
//    return *this;
//}

// python-like function to get the size of a vector
int len( const Vector& v )
{
    return v.mSize;
}

void print( const Vector& v )
{
    assert( v.mSize > 0 );
    for ( int i = 0; i < v.mSize; i++ )
    {
        std::cout << v.mData[i] << '\n';
    }
}
