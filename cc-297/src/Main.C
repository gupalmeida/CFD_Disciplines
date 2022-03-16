#include <iostream>

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

    sol = s.Solve();

    print(sol);

    return 0;
}
