# Adds materials to ErrorEstimationLib
target_sources(ErrorEstimationLib PUBLIC
    TPZHDivErrorEstimateMaterial.cpp
    TPZHDivErrorEstimateMaterial.h
    TPZHDivErrorEstimateDarcyMaterial.cpp
    TPZHDivErrorEstimateDarcyMaterial.h
    TPZHDivErrorEstimateElasticityMaterial.cpp
    TPZHDivErrorEstimateElasticityMaterial.h
    TPZMixedHdivErrorEstimate.cpp
    TPZMixedHdivErrorEstimate
    TPZMatLaplacianHybrid.cpp
    TPZMatLaplacianHybrid.h
    TPZPressureProjection.h
    TPZPressureProjection.cpp
    TPZHybridH1ErrorEstimateMaterial.h
    TPZHybridH1ErrorEstimateMaterial.cpp
    TPZHybridH1HdivFluxRecMaterial.h
    TPZHybridH1HdivFluxRecMaterial.cpp
    TPZHybridH1PressureRecMaterial.h
    TPZHybridH1PressureRecMaterial.cpp
    TPZSteklovMaterial.h
    TPZSteklovMaterial.cpp
)

target_include_directories(ErrorEstimationLib PUBLIC ${CMAKE_CURRENT_LIST_DIR})
