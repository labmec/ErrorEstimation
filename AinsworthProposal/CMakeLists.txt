# This target implements M. Ainsworth proposal using new NeoPZ functionalities.

file(GLOB header *.h)
file(GLOB source *.cpp)

set(ErrorEstimateSRC ${header} ${source})

add_executable(AinsworthProposal ${ErrorEstimateSRC})

target_link_libraries(AinsworthProposal pz)

configure_file(../Meshes/Quad.msh Quad.msh COPYONLY)
configure_file(../Meshes/Cube.msh Cube.msh COPYONLY)
configure_file(../Meshes/BasicMesh.msh BasicMesh.msh COPYONLY)
configure_file(../Meshes/LMesh.msh LMesh.msh COPYONLY)
configure_file(../Meshes/esfera.msh esfera.msh COPYONLY)
