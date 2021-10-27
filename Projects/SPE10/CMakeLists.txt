add_executable(SPE10 mainSPE10.cpp ToolsSPE10.cpp ToolsSPE10.h)
target_link_libraries(SPE10 PRIVATE Tools ErrorEstimationLib)

find_package(libInterpolate REQUIRED)
target_link_libraries(SPE10 PRIVATE libInterpolate::Interpolate)

configure_file(InputData/spe_perm.dat InputData/spe_perm.dat COPYONLY)
