# Set the project name
project (libprimitives)

add_library(
    ${PROJECT_NAME}
    INTERFACE
)

add_library(
    libprimitive
    ALIAS
    ${PROJECT_NAME}
)

target_include_directories(
    ${PROJECT_NAME}
    INTERFACE
    ${PROJECT_SOURCE_DIR}/include
)