# Project to implemented MHM Error Estimation
file(GLOB header *.h)
file(GLOB source *.cpp)

set(ErrorEstimateSRC  ${common_SRC} ${header} ${source})

add_executable(MHMEstimation ${ErrorEstimateSRC})

target_link_libraries(MHMEstimation pz)

configure_file(../Meshes/Quad.msh Quad.msh COPYONLY)
configure_file(../Meshes/Cube.msh Cube.msh COPYONLY)
configure_file(../Meshes/BasicMesh.msh BasicMesh.msh COPYONLY)
configure_file(../Meshes/LMesh.msh LMesh.msh COPYONLY)
configure_file(../Meshes/esfera.msh esfera.msh COPYONLY)
