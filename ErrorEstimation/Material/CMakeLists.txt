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
    TPZEEMatHybridH1ToH1.h
    TPZEEMatHybridH1ToH1.cpp
    TPZEEMatHybridH1ToHDiv.h
    TPZEEMatHybridH1ToHDiv.cpp
	TPZSteklovMaterial.h
	TPZSteklovMaterial.cpp
    )

target_include_directories(ErrorEstimationLib PUBLIC ${CMAKE_CURRENT_LIST_DIR})
