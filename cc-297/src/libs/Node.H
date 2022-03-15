#ifndef NODE_H
#define NODE_H

#include <iostream>
#include <vector>

class Node
{
    private:
        double x_coord;
        double y_coord;
    
    public:
        // Constructors & Destructors
        Node();
        ~Node();

        // Member functions
        void setXCoord( double x );
        void setYCoord( double y );

        double getXCoord( ) const;
        double getYCoord( ) const;
};

#endif
