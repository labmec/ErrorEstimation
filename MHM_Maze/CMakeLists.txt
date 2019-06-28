

set(MHMMAZE mainMaze.cpp ConfigCasesMaze.h meshgen.h meshgen.cpp
)

add_executable(MHM_Maze mainMaze.cpp ${MHMMAZE})
target_link_libraries(MHM_Maze pz ${OpenCV_LIBS})

<<<<<<< HEAD
set(MHMError mainMHMError.cpp ConfigCasesMaze.h meshgen.h meshgen.cpp ${common_SRC})
=======
set(MHMError mainMHMError.cpp ConfigCasesMaze.h meshgen.h meshgen.cpp
        ../ErrorEstimationHDiv/TPZHybridHDivErrorEstimator.h
        ../ErrorEstimationHDiv/TPZHybridHDivErrorEstimator.cpp
        ../ErrorEstimationHDiv/TPZMixedHdivErrorEstimate.h
        ../ErrorEstimationHDiv/TPZMixedHdivErrorEstimate.cpp
        ../ErrorEstimationHDiv/ProblemConfig.h
        ../ErrorEstimationHDiv/TPZHDivErrorEstimateMaterial.h
        ../ErrorEstimationHDiv/TPZHDivErrorEstimateMaterial.cpp
	TPZMHMHDivErrorEstimator.h
	TPZMHMHDivErrorEstimator.cpp
)
>>>>>>> eb2736884061485bc54ace5f8abfe4aab3eb829b

add_executable(MHM_Error ${MHMError})
target_link_libraries(MHM_Error pz ${OpenCV_LIBS})


file(COPY Mazes DESTINATION .)
