#ifndef IOTOOLS_H
#define IOTOOLS_H

#include <sstream>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>

#include "tmpIOTools.H"

class TecWriter
{
public:
    TecWriter();

    std::string fName;

    void SetNumberOfVariables(int nVars);
    void SetTitle(std::string title);
    void SetFileType(std::string fType);
    void SetZoneType(std::string zType);
    void SetDataPacking(std::string zType);
    void SetVariables(std::vector<std::string> variables);
    template<class T> T WriteSolution(T &data);
    void MakeHeader();

private:
    int mNumberOfVariables;
    std::string mTitle, mFileType, mZoneType;
    std::string mDataPacking;
    std::vector<std::string> mVariables;
};

#endif
