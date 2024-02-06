add_library(
    liblinAlg
    src/LinSys.cpp
    src/Matrix.cpp
    src/Vector.cpp
)

target_include_directories(
    liblinAlg
    PUBLIC
    include
)
