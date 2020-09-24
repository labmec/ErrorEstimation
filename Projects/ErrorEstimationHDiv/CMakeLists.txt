# Project to apply error estimation for H(div) simulations
add_executable(HeteroPermeability main_HeteroPermeability.cpp)
target_link_libraries(HeteroPermeability Tools ErrorEstimationLib)

add_executable(HybridHDiv main.cpp)
target_link_libraries(HybridHDiv Tools ErrorEstimationLib)

add_executable(AdaptivityTest mainAdaptivity.cpp)
target_link_libraries(AdaptivityTest Tools ErrorEstimationLib)

add_executable(QuadMesh mainQuadMesh.cpp)
target_link_libraries(QuadMesh Tools ErrorEstimationLib)

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
