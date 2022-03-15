cmake_minimum_required( VERSION 3.15 )

project( cc297 VERSION 1.0 )

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add libraries directory
add_subdirectory( lib )
add_subdirectory( src )

# add the binary tree to the search path for include files
target_include_directories(
    ${PROJECT_NAME} PUBLIC
    "${PROJECT_BINARY_DIR}"
    "${PROJECT_SOURCE_DIR}/src"
    "${PROJECT_SOURCE_DIR}/lib"
)

target_link_libraries( ${PROJECT_NAME} PUBLIC lib )

install(TARGETS ${PROJECT_NAME} DESTINATION bin)
