#include <iostream>
#include <chrono>

#include "Matrix.H"
#include "Vector.H"
// #include "SecondOrderODE.H"
// #include "LinSys.H"
// #include "IOobject.H"

int main()
{
    Matrix m(5,5);
    Vector b(5);

    b(0) = 1.0;
    b(4) = 1.0;

    m.setUpperDiag( 1.0 );
    m.setLowerDiag( 1.0 );
    m.setDiag( -2.0 );
    
    print(m);

    return 0;
}
