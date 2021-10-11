add_library(Methods
        InputTreatment.h
        InputTreatment.cpp
        MeshInit.h
        MeshInit.cpp
        Output.cpp
        Output.h
        Solver.cpp
        Solver.h
        TPZHybMixDiffMaterial.h
        TPZHybMixDiffMaterial.cpp
        ForcingFunction.cpp
        ForcingFunction.h
        )

target_link_libraries(Methods Tools ErrorEstimationLib)
target_include_directories(Methods PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Project to apply error estimation for hybrid H1 approximations
add_executable(HybridH1 main_HybridH1.cpp)
target_link_libraries(HybridH1 Tools ErrorEstimationLib)

# Project to apply error estimation for hybrid H1 approximations
add_executable(HybMixDiff main_HybMixDiff.cpp)
target_link_libraries(HybMixDiff Tools Methods ErrorEstimationLib)

add_executable(HybridH1vsMixed main_HybridH1vsMixed.cpp)
target_link_libraries(HybridH1vsMixed Tools Methods ErrorEstimationLib)

# Project to apply error estimation for hybrid H1 approximations
add_executable(FluxError main_fluxError.cpp)
target_link_libraries(FluxError Tools Methods ErrorEstimationLib)

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
