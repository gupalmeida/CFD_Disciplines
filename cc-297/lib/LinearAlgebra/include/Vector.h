#ifndef VECTOR_H
#define VECTOR_H

#include <cmath>
#include <cassert>
#include <iostream>

class Vector
{
    private:
        double* mData; // data stored in vector
        int mSize; // size of the vector

    public:
        Vector(int size);
        Vector(const Vector& v);
        ~Vector();
        // vector initialization
        void ones();
        // vector info and printing
        int size() const;
        void print() const;
        void print( int i ) const;
        // reading vector element
        double& operator()(int i); // zero-based indexing
        double& operator[](int i); // redundant with operator ()
        // assignment
        Vector& operator=( const Vector& v );
        Vector operator+() const; // unary operator
        Vector operator-() const; // unary operator
        Vector operator+( const Vector& v ) const; // binary operator
        Vector operator-( const Vector& v ) const; // binary operator
        Vector operator*(double a) const;
        //Vector& operator*(double a);

        // Calculate p-norm
        double pNorm( int p=2 ) const;

        friend int len(const Vector& v);
        friend void print(const Vector& v);

};

int len( const Vector& v );
void print(const Vector& v);

#endif
