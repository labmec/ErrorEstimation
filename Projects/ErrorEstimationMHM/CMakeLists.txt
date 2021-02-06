# Project to implemented MHM Error Estimation
add_executable(MHMEstimation main.cpp)
target_link_libraries(MHMEstimation Tools ErrorEstimationLib)

add_executable(MHM_CILAMCE main_CILAMCE.cpp)
target_link_libraries(MHM_CILAMCE Tools ErrorEstimationLib)

add_executable(MethodsExamples mainExampleMethods.cpp)
target_link_libraries(MethodsExamples Tools ErrorEstimationLib)

configure_file(../Meshes/Quad.msh Quad.msh COPYONLY)
configure_file(../Meshes/Cube.msh Cube.msh COPYONLY)
configure_file(../Meshes/BasicMesh.msh BasicMesh.msh COPYONLY)
configure_file(../Meshes/LMesh.msh LMesh.msh COPYONLY)
configure_file(../Meshes/esfera.msh esfera.msh COPYONLY)
configure_file(../Meshes/Circular.msh Circular.msh COPYONLY)
