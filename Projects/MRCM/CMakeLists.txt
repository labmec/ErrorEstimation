# Project to apply error estimation for hybrid H1 approximations
add_executable(MRCM main_MRCM.cpp)
target_link_libraries(MRCM Tools ErrorEstimationLib)
