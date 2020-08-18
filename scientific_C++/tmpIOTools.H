#ifndef TEMPIOTOOLS_H
#define TEMPIOTOOLS_H

#include <sstream>
#include <string>
#include <iostream>
#include <fstream>
#include <cassert>


//template<class T> T GetData(std::string varName, std::string fileName);]
template<class T> T GetData(std::string varName, std::string fileName)
{
    T varValue;
    std::ifstream inputFile(fileName);
    assert(inputFile.is_open());
    std::string line;
    
    while ( getline(inputFile,line) )
    {
        if ( ("@" + varName) == line )
        {
            break;
        }
    }

    inputFile >> line;

    std::stringstream value(line);

    value >> varValue;

    inputFile.close();

    return varValue;
}

#endif
