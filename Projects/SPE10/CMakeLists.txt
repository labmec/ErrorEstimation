add_executable(SPE10 mainSPE10.cpp)
target_link_libraries(SPE10 PRIVATE Tools ErrorEstimationLib)

configure_file(InputData/spe_perm.dat InputData/spe_perm.dat COPYONLY)
