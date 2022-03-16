#ifndef MESH_H
#define MESH_H

#include <iostream>
#include <cmath>
#include <vector>
#include <cassert>

#include "Macros.H"
#include "Node.H"

class Mesh
{
    private:
        std::vector<Node> mNodes;

    public:
        Mesh( double xMin, double xMax, int numNodes );
};

#endif
