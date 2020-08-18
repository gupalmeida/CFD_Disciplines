#ifndef VECTOR_H
#define VECTOR_H

#include <iostream>
#include <vector>
#include <cmath>

class Vector
{
private:
    double* mData;
    int mSize;

public:
    Vector(const Vector& otherVector);
    Vector(int size);
    ~Vector();
    int GetSize() const;
    double& operator[] (int i);
    double& operator() (int i);
    double Read(int i) const;
    Vector& operator=(const Vector& otherVector);
    Vector operator+() const;
    Vector operator-() const;
    Vector operator+(const Vector& v1) const;
    Vector operator-(const Vector& v1) const;
    // scalar product
    Vector operator* (double a) const;

    // p-norm
    double ComputePNorm(int p=2) const;
    friend int len(const Vector& v);
    friend void print(const Vector& v);

};

int len(const Vector& v);
void print(const Vector& v);

#endif
