#ifndef MESH_H
#define MESH_H

#include <iostream>
#include <cmath>
#include <vector>
#include <cassert>
using namespace std;

#include "Node.H"

class Mesh
{
    public:
        friend class BvpOde;

    private:
        vector<Node> mNodes;

    public:
        Mesh( double xMin, double xMax, int numNodes );
};

#endif
