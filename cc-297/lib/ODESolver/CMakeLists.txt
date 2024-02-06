add_library(
    libode
    src/BoundaryConditions.cpp
    src/BvpOde.cpp
    src/Mesh.cpp
)

target_include_directories(
    libode
    PUBLIC
    include
)

target_link_libraries(
    libode
    liblinAlg
    libioobj
)
