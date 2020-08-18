#include <cassert>
#include "IOTools.H"

// TecWrite class for writing output file in Tecplot ASCII data format
TecWriter::TecWriter()
{
    mTitle = "SIMULATION";
    mFileType = "FULL";
    mZoneType = "ORDERED";
    mDataPacking = "POINT";
    fName = "output.dat";
}

void TecWriter::SetNumberOfVariables(int nVars)
{
    assert(nVars > 0);
    mNumberOfVariables = nVars;
}

void TecWriter::SetTitle(std::string title)
{
    mTitle = title;
}

void TecWriter::SetFileType(std::string fType)
{
    mFileType = fType;
}

void TecWriter::SetZoneType(std::string zType)
{
    mZoneType = zType;
}

void TecWriter::SetDataPacking(std::string dPacking)
{
    mDataPacking = dPacking;
}

void TecWriter::SetVariables(std::vector<std::string> variables)
{
    assert( variables.size() != 0 );
    mVariables = variables;
}

void TecWriter::MakeHeader()
{
    std::ofstream outputFile(fName);
    assert( outputFile.is_open() );
    int nx, ny, nz;

    // getting number of points in each direction
    nx = GetData<int>("nx","input.inp");
    ny = GetData<int>("ny","input.inp");
    nz = GetData<int>("nz","input.inp");

    // write file title
    outputFile << "TITLE = ";
    outputFile << "\"" + mTitle + "\"" << "\n";

    // write file type
    outputFile << "FILETYPE = ";
    outputFile << "\"" + mFileType + "\"" << "\n";

    // write variables
    outputFile << "VARIABLES = ";
    for ( int i = 0; i < mVariables.size(); i++)
    {
        outputFile << mVariables[i] + ",";
    }
    outputFile << "\n";

    // write zone data
    outputFile << "ZONE\n";
    outputFile << "ZONETYPE = ";
    outputFile << "\"" + mZoneType + "\"" << "\n";
    outputFile << "I = " << nx << ", ";
    outputFile << "J = " << ny << ", ";
    outputFile << "K = " << nz << "\n";

    outputFile <<"ZONETYPE = " << mZoneType << ", ";
    outputFile <<"DATAPACKING = " << mDataPacking << "\n";
    
    // close file
    outputFile.close();
}

template<class T> T TecWriter::WriteSolution(T &data)
{
    std::ofstream outputFile(fName,std::ios::app);
    assert( outputFile.is_open() );
    outputFile.close();
}

// Explicit instantiation of template functions
template int TecWriter::WriteSolution(int &data);
template double TecWriter::WriteSolution(double &data);
template std::string TecWriter::WriteSolution(std::string &data);
