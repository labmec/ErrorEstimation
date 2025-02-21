# This project estimates the error of a flow around a NACA profile
set(ErrorNACASRC
    ErrorNaca.cpp tpznacaprofile.h tpznacaprofile.cpp
)

add_executable(ErrorNACA ${ErrorNACASRC})
target_link_libraries(ErrorNACA NeoPZ::pz)

add_executable(2DSqrt 2DSqrtProblem.cpp)
target_link_libraries(2DSqrt NeoPZ::pz)

add_executable(Lshaped LshapedProblem.cpp)
target_link_libraries(Lshaped NeoPZ::pz)

add_executable(Constrastingpermeability ConstrastingpermeabilityProblem.cpp)
target_link_libraries(Constrastingpermeability NeoPZ::pz)

# copy a file from source directory to binary directory
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/naca.msh ${CMAKE_CURRENT_BINARY_DIR}/naca.msh COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/2DSqrtProblem.msh ${CMAKE_CURRENT_BINARY_DIR}/2DSqrtProblem.msh COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/LshapedProblem.msh ${CMAKE_CURRENT_BINARY_DIR}/LshapedProblem.msh COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/ConstrastingpermeabilityProblem.msh ${CMAKE_CURRENT_BINARY_DIR}/ConstrastingpermeabilityProblem.msh COPYONLY)