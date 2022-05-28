# Set the project name
project (liblinearAlgebra)

# Add a library with the above sources
add_library(
    ${PROJECT_NAME}
    src/LinSys.C
    src/Matrix.C
    src/Vector.C
)

add_library(
    liblinAlg
    ALIAS
    ${PROJECT_NAME}
)

target_include_directories(
    ${PROJECT_NAME}
    PUBLIC
    ${PROJECT_SOURCE_DIR}/include
)