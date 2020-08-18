#include <cassert>
#include "Vector.H"

Vector::Vector(const Vector& otherVector)
{
    mSize = otherVector.GetSize();
    mData = new double [mSize];
    for (int i = 0; i < mSize; i++)
    {
        mData[i] = otherVector.mData[i];
    }
}

Vector::Vector(int size)
{
    assert( size > 0 );
    mSize = size;
    mData = new double [mSize];
    for (int i= 0; i < mSize; i++)
    {
        mData[i] = 0.0;
    }
}

Vector::~Vector()
{
    // free memory of object instantiation
    //std::cout << "Deleting Vector object\n";

    delete[] mData;
}

int Vector::GetSize() const
{
    return mSize;
}

double& Vector::operator[] (int i)
{
    assert( i > -1);
    assert( i < mSize );
    return mData[i];
}

double& Vector::operator() (int i)
{
    assert( i > -1);
    assert( i < mSize );
    return mData[i];
}

double Vector::Read(int i) const
{
    return mData[i];
}

Vector& Vector::operator= (const Vector& otherVector)
{
    assert( mSize == otherVector.mSize );
    for (int i = 0; i < mSize; i++)
    {
        mData[i] = otherVector.mData[i];
    }

    return *this;
}

Vector Vector::operator+ () const
{
    Vector v(mSize);

    for (int i = 0; i < mSize; i++)
    {
        v[i] = mData[i];
    }

    return v;
}

Vector Vector::operator- () const
{
    Vector v(mSize);

    for (int i = 0; i < mSize; i++)
    {
        v[i] = -mData[i];
    }

    return v;
}

Vector Vector::operator+ (const Vector& v1) const
{
    assert( mSize == v1.mSize );
    Vector v(mSize);

    for (int i = 0; i < mSize; i++)
    {
        v[i] = mData[i] + v1.mData[i];
    }

    return v;
}

Vector Vector::operator- (const Vector& v1) const
{
    assert( mSize == v1.mSize );
    Vector v(mSize);

    for (int i = 0; i < mSize; i++)
    {
        v[i] = mData[i] - v1.mData[i];
    }

    return v;
}

Vector Vector::operator* (double a) const
{
    Vector v(mSize);

    for (int i = 0; i < mSize; i++)
    {
        v[i] = a * mData[i];
    }

    return v;
}

double Vector::ComputePNorm(int p) const
{
    double normVal, sum = 0.0;

    for (int i = 0; i < mSize; i++)
    {
        sum += pow( fabs(mData[i]), p );
    }

    normVal = pow( sum, (1.0 / ((double) p)) );

    return normVal;
}

int len(const Vector& v)
{
    return v.mSize;
}

void print(const Vector& v)
{
    std::cout << "[ ";

    for (int i = 0; i < v.mSize; i++)
    {
        std::cout << v.mData[i] << " ";
    }

    std::cout << "]\n";
}
