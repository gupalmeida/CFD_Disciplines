#ifndef IOOBJECT_H
#define IOOBJECT_H

#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>
#include <cassert>

#include "Aux.H"
#include "Vector.H"
#include "Matrix.H"

using namespace std;

class IOobject
{
    public:
        void readSetup( );
};

#endif
