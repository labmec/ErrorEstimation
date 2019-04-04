# This project works with error estimation on H1 and H(div) spaces.

set(ErrorEstimateSRCGmesh main.cpp TPZMixedErrorEstimate.h TPZMixedErrorEstimate.cpp TPZPostProcessError.h
        TPZPostProcessError.cpp ProblemConfig.h)

add_executable(ErrorEstimateH1GMsh ${ErrorEstimateSRCGmesh})
target_link_libraries(ErrorEstimateH1GMsh pz)

configure_file(Quad.msh Quad.msh COPYONLY)
configure_file(Cube.msh Cube.msh COPYONLY)
