#ifndef FDGRID_H
#define FDGRID_H

#include <cmath>
#include <vector>
#include "Node.H"

class FDGrid
{
public:
    friend class BvpODE;

private:
    std::vector<Node> mNodes;

public:
    FDGrid(int numNodes, double xMin, double xMax);
};

#endif
