# Adds materials to ErrorEstimationLib
target_sources(ErrorEstimationLib PUBLIC
    TPZHDivErrorEstimateMaterial.cpp
    TPZHDivErrorEstimateMaterial.h
    TPZMHMHDivErrorEstimatorMaterial.cpp
    TPZMHMHDivErrorEstimatorMaterial.h
    TPZMatLaplacianHybrid.cpp
    TPZMatLaplacianHybrid.h
)

target_include_directories(ErrorEstimationLib PUBLIC ${CMAKE_CURRENT_LIST_DIR})
