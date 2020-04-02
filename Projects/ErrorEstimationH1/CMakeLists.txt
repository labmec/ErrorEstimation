# This project works with error estimation on H1 and H(div) spaces.
set(ErrorEstimateSRCGmesh
    main.cpp
    TPZMixedErrorEstimate.h
    TPZMixedErrorEstimate.cpp
    TPZPostProcessError.h
    TPZPostProcessError.cpp
    ProblemConfig.h
)

add_executable(ErrorEstimateH1GMsh ${ErrorEstimateSRCGmesh})
target_link_libraries(ErrorEstimateH1GMsh Tools ErrorEstimationLib)

configure_file(../Meshes/Quad.msh Quad.msh COPYONLY)
configure_file(../Meshes/Cube.msh Cube.msh COPYONLY)
configure_file(../Meshes/BasicMesh.msh BasicMesh.msh COPYONLY)
configure_file(../Meshes/LMesh.msh LMesh.msh COPYONLY)
configure_file(../Meshes/esfera.msh esfera.msh COPYONLY)
