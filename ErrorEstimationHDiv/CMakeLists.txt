# Project to apply hp-adaptive process to Poisson3D material on hexahedra with a shock problem

file(GLOB header *.h)
file(GLOB source *.cpp)

set(ErrorEstimateSRC ${header} ${source})

add_executable(ErrorEstimateHDiv ${ErrorEstimateSRC})

target_link_libraries(ErrorEstimateHDiv pz)

configure_file(BasicMesh.msh BasicMesh.msh COPYONLY)
