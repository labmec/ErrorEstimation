# This project estimates the error of a flow around a NACA profile
set(ErrorNACASRC
    ErrorNacaH1HDiv.cpp tpznacaprofile.h tpznacaprofile.cpp
)

add_executable(ErrorNACAH1HDiv ${ErrorNACASRC})
target_link_libraries(ErrorNACAH1HDiv NeoPZ::pz)

# copy a file from source directory to binary directory
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/naca.msh ${CMAKE_CURRENT_BINARY_DIR}/naca.msh COPYONLY)