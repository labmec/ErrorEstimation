# Project to apply error estimation for Hdiv reconstrcuction
set(
        common_SRC
        ProblemConfig.h
        Tools.h
        Tools.cpp
        TPZHybridHDivErrorEstimator.h
        TPZHybridHDivErrorEstimator.cpp
        TPZMixedHdivErrorEstimate.h
        TPZMixedHdivErrorEstimate.cpp
        TPZHDivErrorEstimatorH1.h
        TPZHDivErrorEstimatorH1.cpp
        TPZHDivErrorEstimateMaterial.h
        TPZHDivErrorEstimateMaterial.cpp
)

add_executable(HeteroPermeability ${common_SRC} main_HeteroPermeability.cpp)

target_link_libraries(HeteroPermeability pz)


add_executable(HybridMain ${common_SRC} New_main.cpp)

target_link_libraries(HybridMain pz)


add_executable(NoHybridMain ${common_SRC} main.cpp)

target_link_libraries(NoHybridMain pz)

add_executable(H1Reconstruction ${common_SRC} main_H1.cpp)

target_link_libraries(H1Reconstruction pz)


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
