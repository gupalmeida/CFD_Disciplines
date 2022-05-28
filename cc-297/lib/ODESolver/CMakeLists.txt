# Set the project name
project (libodeSolver)

# Add a library with the above sources
add_library(
    ${PROJECT_NAME}
    src/BoundaryConditions.C
    src/BvpOde.C
    src/Mesh.C
)

add_library(
    libode
    ALIAS
    ${PROJECT_NAME}
)

target_include_directories(
    ${PROJECT_NAME}
    PUBLIC
    ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(
    ${PROJECT_NAME}
    liblinAlg
    libioobj
)