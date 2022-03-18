#ifndef MESH_H
#define MESH_H

#include <iostream>
#include <cmath>
#include <vector>
#include <cassert>

#include "Aux.H"
#include "Node.H"

class Mesh
{
    private:
        vector<Node> mNodes;

    public:
        Mesh( double xMin, double xMax, int numNodes );
};

#endif
