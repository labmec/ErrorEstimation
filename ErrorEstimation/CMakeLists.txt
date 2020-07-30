# This file creates a library with error estimation files
add_library(ErrorEstimationLib
    # Error estimation
    ProblemConfig.h
    TPZHDivErrorEstimatorH1.cpp
    TPZHDivErrorEstimatorH1.h
    TPZHybridHDivErrorEstimator.cpp
    TPZHybridHDivErrorEstimator.h
    TPZHybridH1ErrorEstimator.cpp
    TPZHybridH1ErrorEstimator.h
    TPZMHMHDivErrorEstimator.cpp
    TPZMHMHDivErrorEstimator.h
    TPZMixedHdivErrorEstimate.cpp
    TPZMixedHdivErrorEstimate.h
    # Related tools
    TPZCreateMultiphysicsSpace.cpp
    TPZCreateMultiphysicsSpace.h
    TPZMFSolutionTransfer.cpp
    TPZMFSolutionTransfer.h

)
target_link_libraries(ErrorEstimationLib pz)
target_include_directories(ErrorEstimationLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Add material files to ErrorEstimationLib
add_subdirectory(Material)
