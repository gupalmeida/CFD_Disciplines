project( cc297 )

# Create the executable
add_executable(
    ${PROJECT_NAME}
    # Main.C
    Driver.C
)

# Link the static library from subproject1 using its alias sub::lib1
# Link the header only library from subproject2 using its alias sub::lib2
# This will cause the include directories for that target to be added to this project
target_link_libraries(
    ${PROJECT_NAME}
    liblinAlg
    libode
    libprimitive
)