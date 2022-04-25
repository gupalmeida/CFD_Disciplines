cmake_minimum_required( VERSION 3.15 )

project( cc297 VERSION 1.0 )

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_BUILD_TYPE Debug)

# Add libraries directory
add_subdirectory( lib )
add_subdirectory( src )

# Set root directory of source code
set(
    BUILD_SRC_DIR
    "/home/vostrocfd/git/CFD_Disciplines/cc-297"
    CACHE PATH
    "Sets root directory of source code."
    FORCE
)

# add the binary tree to the search path for include files
target_include_directories(
    ${PROJECT_NAME} PUBLIC
    "${PROJECT_BINARY_DIR}"
    "${BUILD_SRC_DIR}/src"
    "${BUILD_SRC_DIR}/lib"
)

target_link_libraries( ${PROJECT_NAME} PUBLIC lib )

install(TARGETS ${PROJECT_NAME} DESTINATION bin)
