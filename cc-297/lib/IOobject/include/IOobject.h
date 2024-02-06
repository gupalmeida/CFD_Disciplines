#ifndef IOOBJECT_H
#define IOOBJECT_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <cassert>
#include <string>

namespace fs = std::filesystem;

void fileExists(const std::string& path);

template<typename T> T readInput( std::string varname, std::string& infile );

#endif
