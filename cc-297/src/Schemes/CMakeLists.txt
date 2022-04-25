set(
    FILES
    "ConvectiveSchemes.H"
    "ConvectiveSchemes.C"
)

target_sources(
    "${PROJECT_NAME}"
    PRIVATE
    ${FILES}
)
