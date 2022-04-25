set(
    SUBDIRECTORIES
    "Schemes"
)

set(
    FILES
    "Driver.C"
    "IOobject.H"
    "Mesh.C"
    "Mesh.H"
    "BoundaryConditions.H"
    "BoundaryConditions.C"
    "BvpOde.H"
    "BvpOde.C"
    "SecondOrderODE.H"
)

#target_sources(
#    "${PROJECT_NAME}"
#    PRIVATE
#    ${FILES}
#)

# Add all subdirectories
foreach( VAR ${SUBDIRECTORIES} )
    add_subdirectory( ${VAR} )
endforeach()
    
target_sources(
    "${PROJECT_NAME}"
    PRIVATE
    ${FILES}
)
install(FILES ${FILES} DESTINATION "${BUILD_SRC_DIR}/build")
#add_executable(
#    ${PROJECT_NAME}
#    ${FILES}
#)

#message( STATUS "TARGETS: ${FILES}")
