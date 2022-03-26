#ifndef IOOBJECT_H
#define IOOBJECT_H

#include <iostream>
#include <fstream>
#include <cassert>
#include <string>

#include "Aux.H"

using namespace std;

template<class T>
T readInput( string varname, string filename )
{
    // declaring variable to get as input from varname
    T var;

    ifstream readFile;
    readFile.open( filename );

    string meof;
    bool found = 0;
    while ( readFile >> meof && !found )
    {
        if( meof == varname )
        {
            found = 1;
            readFile.ignore();
            readFile >> var;
        }
    }

    readFile.close();

    return var;
}

#endif
