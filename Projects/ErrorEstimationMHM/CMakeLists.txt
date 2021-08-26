# Project to implemented MHM Error Estimation
add_executable(MHM_Estimation main_MHM.cpp)
target_link_libraries(MHM_Estimation Tools ErrorEstimationLib)

add_executable(MHM_Journal main_Journal.cpp)
target_link_libraries(MHM_Journal Tools ErrorEstimationLib)

add_executable(MHM_RefSol main_RefSolValidation.cpp)
target_link_libraries(MHM_RefSol Tools ErrorEstimationLib)
