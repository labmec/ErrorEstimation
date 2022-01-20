add_executable(SPE10Adaptivity mainSPE10.cpp ToolsSPE10.cpp ToolsSPE10.h TPZMultiscaleGridGen2D.cpp TPZMultiscaleGridGen2D.h)
target_link_libraries(SPE10Adaptivity PRIVATE Tools ErrorEstimationLib)

add_executable(SPE10RefSol mainSPE10RefSol.cpp ToolsSPE10.cpp ToolsSPE10.h TPZMultiscaleGridGen2D.cpp TPZMultiscaleGridGen2D.h)
target_link_libraries(SPE10RefSol PRIVATE Tools ErrorEstimationLib)

configure_file(InputData/spe_perm.dat InputData/spe_perm.dat COPYONLY)
