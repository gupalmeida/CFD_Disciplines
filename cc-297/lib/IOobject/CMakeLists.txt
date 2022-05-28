# Set the project name
project (libioobjects)

add_library(
    ${PROJECT_NAME}
    INTERFACE
)

add_library(
    libioobj
    ALIAS
    ${PROJECT_NAME}
)

target_include_directories(
    ${PROJECT_NAME}
    INTERFACE
    ${PROJECT_SOURCE_DIR}/include
)