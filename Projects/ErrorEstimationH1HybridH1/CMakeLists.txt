# This project works with error estimation on H1 and H(div) spaces.
set(ErrorEstimateSRC
    main.cpp
    TPZMixedErrorEstimate2.h
    TPZMixedErrorEstimate2.cpp
    TPZPostProcessError2.h
    TPZPostProcessError2.cpp
    ProblemConfig.h
#    ProblemConfig.cpp
)

add_executable(ErrorEstimateH1HybridH1 ${ErrorEstimateSRC})
target_link_libraries(ErrorEstimateH1HybridH1 Tools Methods ErrorEstimationLib)

configure_file(../Meshes/Quad.msh Quad.msh COPYONLY)
configure_file(../Meshes/Cube.msh Cube.msh COPYONLY)
configure_file(../Meshes/BasicMesh.msh BasicMesh.msh COPYONLY)
configure_file(../Meshes/LMesh.msh LMesh.msh COPYONLY)
configure_file(../Meshes/esfera.msh esfera.msh COPYONLY)
