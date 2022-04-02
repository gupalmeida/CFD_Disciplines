#include <iostream>
#include <chrono>

#include "Matrix.H"
#include "Vector.H"
#include "SecondOrderODE.H"
#include "LinSys.H"
#include "IOobject.H"

int main()
{
    Matrix m(5,5);
    Vector b(5);

    b(0) = 1.0;
    b(4) = 1.0;

    m.setUpperDiag( 1.0 );
    m.setLowerDiag( 1.0 );
    m.setDiag( -2.0 );

    LinSys s(m,b);
    Vector sol(5);

    // measuring time
    auto start = chrono::steady_clock::now();

    sol = s.Solve();

    print( sol );

    auto end = chrono::steady_clock::now();
    cout
        << "Elapsed Solver Time [s]: "
        << chrono::duration_cast<chrono::nanoseconds>(end - start).count()
        << "\n";

    double imax;
    imax = readInput<double>( "JMAX", SETUP_FILE );
    cout << "Value from input file: " << imax << "\n";

    return 0;
}
