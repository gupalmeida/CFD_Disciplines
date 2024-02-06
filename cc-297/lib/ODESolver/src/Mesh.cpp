#include "Mesh.h"

Mesh::Mesh(double xMin, double xMax, int numNodes)
{
    double dx = (xMax - xMin)/(double)(numNodes-1);

    for (int i = 0; i < numNodes; i++)
    {
        Node node;
        node.xCoord = xMin + (double)(i * dx);
        mNodes.push_back(node);
    }

    assert( mNodes.size() == numNodes );
}
