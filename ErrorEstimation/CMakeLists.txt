# This file creates a library with error estimation files
add_library(ErrorEstimationLib
    # Error estimation
    ProblemConfig.h
    ProblemConfig.cpp
    TPZHDivErrorEstimator.cpp
    TPZHDivErrorEstimator.h
    TPZHybridH1ErrorEstimator.cpp
    TPZHybridH1ErrorEstimator.h
    TPZHybridH1CreateHDivReconstruction.h
    TPZHybridH1CreateHDivReconstruction.cpp
    TPZHybridH1CreateH1Reconstruction.h
    TPZHybridH1CreateH1Reconstruction.cpp
    TPZHybridH1ReconstructionBase.h
    TPZHybridH1ReconstructionBase.cpp
    TPZMHMHDivErrorEstimator.cpp
    TPZMHMHDivErrorEstimator.h
	TPZElasticityErrorEstimator.h
	TPZElasticityErrorEstimator.cpp
    # Related tools
    TPZCreateMultiphysicsSpace.cpp
    TPZCreateMultiphysicsSpace.h
    TPZMFSolutionTransfer.cpp
    TPZMFSolutionTransfer.h
    )

target_link_libraries(ErrorEstimationLib PRIVATE NeoPZ::pz)
target_include_directories(ErrorEstimationLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${PZ_INCLUDE_DIRS})

# Add material files to ErrorEstimationLib
add_subdirectory(Material)
