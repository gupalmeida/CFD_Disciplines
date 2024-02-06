#include "IOobject.h"

void fileExists(const std::string& path){
  if( fs::exists(path) ){
    std::cout << "Input file found.\n";
  }
  else{
    std::cout << "Input file not found.\n";
  }
}

template<typename T> T readInput( std::string varname, std::string& infile ){
  // reads variable value from input file
  fileExists( infile );
  // declaring variable to get as input from varname
  T var;

  std::ifstream readFile( infile );
  std::string line;
  bool found = false;
  while ( std::getline(readFile,line) ){
    if( line == varname ){
      found = true;
      break;
    }
  }

  if( found && std::getline(readFile,line)){
    std::istringstream iss(line);
    if( !(iss >> var) ){
      std::cerr << "Failed to convert value to desired value.\n";
      throw std::runtime_error( "Conversion failed for variable: "+varname );
    }
  }
  else if( found ){
    std::cout << "Variable found in input file but no value assigned.\n";
  }
  else{
    std::cout << "Variable not found in input file.\n";
  }

  readFile.close();

  return var;
}

  // void writeSolutionToFile( Vector& v, string fname )
  // {
  //    ofstream outputFile( fname );

  // //    int size = len( v );
  //    for ( int i=0; i < v.size(); i++ )
  //    {
  //        outputFile<< v(i) << "\n";
  //    }

  //    outputFile.close();
  // }
//
