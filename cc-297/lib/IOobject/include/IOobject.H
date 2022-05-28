#ifndef IOOBJECT_H
#define IOOBJECT_H

#include <iostream>
#include <fstream>
#include <cassert>
#include <string>
using namespace std;

template<class T>
T readInput( string varname, string filename )
{
    // reads variable value from input file

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

//void writeSolutionToFile( Vector& v, string fname )
//{
//    ofstream outputFile( fname );
//
//    int size = len( v );
//    for ( int i=0; i < v.size(); i++ )
//    {
//        outputFile<< v(i) << "\n";
//    }
//
//    outputFile.close();
//}

#endif
