#include <iostream>
#include <chrono>

#include "Matrix.h"
#include "Vector.h"
// #include "SecondOrderODE.h"
// #include "LinSys.h"
// #include "IOobject.h"

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
