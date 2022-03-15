#include <iostream>

#include "Node.H"

int main()
{
    Node n;

    n.setXCoord( 12.0 );

    std::cout << "Value of x is: " << n.getXCoord() << "\n";

    return 0;
}
