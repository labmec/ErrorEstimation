# Project to apply error estimation for hybrid H1 approximations
add_executable(HybridH1 main_HybridH1.cpp)
target_link_libraries(HybridH1 Tools ErrorEstimationLib)

# Project to apply error estimation for hybrid H1 approximations
add_executable(HybridH1vsMixed main_HybridH1vsMixed.cpp)
target_link_libraries(HybridH1vsMixed Tools ErrorEstimationLib)

configure_file(../Meshes/Quad.msh Quad.msh COPYONLY)
configure_file(../Meshes/Cube.msh Cube.msh COPYONLY)
configure_file(../Meshes/BasicMesh.msh BasicMesh.msh COPYONLY)
configure_file(../Meshes/LMesh.msh LMesh.msh COPYONLY)
configure_file(../Meshes/esfera.msh esfera.msh COPYONLY)
configure_file(../Meshes/Circular.msh Circular.msh COPYONLY)
configure_file(../Meshes/LCircle.msh LCircle.msh COPYONLY)
configure_file(../Meshes/MeshHetero.msh MeshHetero.msh COPYONLY)
configure_file(../Meshes/MalhaSteklov.msh MalhaSteklov.msh COPYONLY)
configure_file(../Meshes/Galvismesh.msh Galvismesh.msh COPYONLY)
configure_file(../Meshes/esfera2.msh Galvismesh.msh COPYONLY)
