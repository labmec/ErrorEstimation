# This project estimates the error of a flow around a NACA profile
set(ErrorSRC
    ErrorH1HDiv.cpp tpznacaprofile.h tpznacaprofile.cpp
)

add_executable(ErrorH1HDiv ${ErrorSRC})
target_link_libraries(ErrorH1HDiv NeoPZ::pz)

# copy a file from source directory to binary directory
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/naca.msh ${CMAKE_CURRENT_BINARY_DIR}/naca.msh COPYONLY)