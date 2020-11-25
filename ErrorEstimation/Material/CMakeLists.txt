# Adds materials to ErrorEstimationLib
target_sources(ErrorEstimationLib PUBLIC
    TPZHDivErrorEstimateMaterial.cpp
    TPZHDivErrorEstimateMaterial.h
    TPZMHMHDivErrorEstimatorMaterial.cpp
    TPZMHMHDivErrorEstimatorMaterial.h
    TPZMatLaplacianHybrid.cpp
    TPZMatLaplacianHybrid.h
    TPZPressureProjection.h
    TPZPressureProjection.cpp
        TPZHybridH1ErrorEstimateMaterial.h
        TPZHybridH1ErrorEstimateMaterial.cpp
    )

target_include_directories(ErrorEstimationLib PUBLIC ${CMAKE_CURRENT_LIST_DIR})
