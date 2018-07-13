# Project to apply hp-adaptive process to Poisson3D material on hexahedra with a shock problem

set(ErrorEstimateSRC main.cpp problem.h Solution2.cpp TPZCreateMesh.h TPZCreateMesh.cpp TPZMixedErrorEstimate.h TPZMixedErrorEstimate.cpp TPZPostProcessError.h TPZPostProcessError.cpp)

set(ErrorEstimateSRCGmesh mainGMesh.cpp TPZMixedErrorEstimate.h TPZMixedErrorEstimate.cpp TPZPostProcessError.h TPZPostProcessError.cpp)


add_executable(ErrorEstimateH1 ${ErrorEstimateSRC})

target_link_libraries(ErrorEstimateH1 pz)

add_executable(ErrorEstimateH1GMsh ${ErrorEstimateSRCGmesh})

target_link_libraries(ErrorEstimateH1GMsh pz)

configure_file(Quad.geo Quad.geo COPYONLY)
configure_file(Cube.geo Cube.geo COPYONLY)
