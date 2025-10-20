# This project works with error estimation on H1 and H(div) spaces.
set(ErrorEstimateSRC
    main_didactic.cpp
    TPZMixedErrorEstimate2.h
    TPZMixedErrorEstimate2.cpp
    TPZPostProcessError2.h
    TPZPostProcessError2.cpp
    ProblemConfig.h
#    ProblemConfig.cpp
)

set(HybridSBemSRC
    main_hybridSBFem.cpp
    TPZBuildSBFemHybrid.h
    TPZBuildSBFemHybrid.cpp
    TPZPostProcessErrorSBFem.h
    TPZPostProcessErrorSBFem.cpp
    TPZSBFemElementGroupPostProcess.h
    TPZSBFemElementGroupPostProcess.cpp
    TPZGeoPatch.h
    TPZGeoPatch.cpp
    )

add_executable(ErrorEstimateSimpleH1HybridH1 ${ErrorEstimateSRC})
target_link_libraries(ErrorEstimateSimpleH1HybridH1 Tools Methods ErrorEstimationLib)

find_package(nlohmann_json 3.2.0 REQUIRED)

add_executable(HybridSBem ${HybridSBemSRC})
target_link_libraries(HybridSBem PRIVATE Tools Methods ErrorEstimationLib nlohmann_json::nlohmann_json)

configure_file(hybrid.json hybrid.json COPYONLY)
configure_file(SqrtSingular.json SqrtSingular.json COPYONLY)
configure_file(../Meshes/Quad.msh Quad.msh COPYONLY)
configure_file(../Meshes/QuadSqrt.msh QuadSqrt.msh COPYONLY)
configure_file(../Meshes/Cube.msh Cube.msh COPYONLY)
configure_file(../Meshes/BasicMesh.msh BasicMesh.msh COPYONLY)
configure_file(../Meshes/LMesh.msh LMesh.msh COPYONLY)
configure_file(../Meshes/esfera.msh esfera.msh COPYONLY)
