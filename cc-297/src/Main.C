#include <iostream>
#include <chrono>

#include "LinSys.H"

int main()
{
    Matrix m(5,5);
    Vector b(5);

    b.ones();
    b(0) = 3.0;
    b(4) = 3.0;

    m.ones();
    m.setUpperDiag( 2.0 );
    m.setLowerDiag( 3.0 );
    m.setDiag( 5.0 );

    print( b );
    std::cout << "\n\n";
    print( m );

    LinSys s(m,b);
    Vector sol(5);

    // measuring time
    auto start = std::chrono::steady_clock::now();

    sol = s.Solve();

    auto end = std::chrono::steady_clock::now();
    std::cout
        << "Elapsed Solver Time [s]: "
        << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()
        << "\n";

    print(sol);

    return 0;
}
