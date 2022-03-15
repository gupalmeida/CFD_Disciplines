#include "Node.H"

// Constructor
Node::Node()
{
    x_coord = -999.0;
    y_coord = -999.0;
}

void Node::setXCoord( double x )
{
    x_coord = x;
}

void Node::setYCoord( double y )
{
    y_coord = y;
}

double Node::getXCoord( ) const
{
    return x_coord;
}

double Node::getYCoord( ) const
{
    return y_coord;
}

// Destructor
Node::~Node()
{
}
