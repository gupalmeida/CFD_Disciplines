cmake_minimum_required( VERSION 3.10 )

project( cc297 VERSION 1.0 )

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Add libraries directory
add_subdirectory( libs )

# Add the executable
add_executable( ${PROJECT_NAME} main.C )

target_link_libraries( ${PROJECT_NAME} PUBLIC libs )

# add the binary tree to the search path for include files
target_include_directories(
    ${PROJECT_NAME} PUBLIC
    "${PROJECT_BINARY_DIR}"
    "${PROJECT_SOURCE_DIR}/libs"
)
