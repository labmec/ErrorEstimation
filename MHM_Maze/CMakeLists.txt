

set(MHMMAZE mainMaze.cpp ConfigCasesMaze.h meshgen.h meshgen.cpp
)

add_executable(MHM_Maze mainMaze.cpp ${MHMMAZE})
target_link_libraries(MHM_Maze pz ${OpenCV_LIBS})

set(MHMError mainMHMError.cpp ConfigCasesMaze.h meshgen.h meshgen.cpp
        ../TPZHybridHDivErrorEstimator.h
        ../TPZHybridHDivErrorEstimator.cpp
        ../TPZMixedHdivErrorEstimate.h
        ../TPZMixedHdivErrorEstimate.cpp
        ../ProblemConfig.h
        ../TPZHDivErrorEstimateMaterial.h
        ../TPZHDivErrorEstimateMaterial.cpp
	../TPZMHMHDivErrorEstimator.h
	../TPZMHMHDivErrorEstimator.cpp
)

add_executable(MHM_Error ${MHMError})
target_link_libraries(MHM_Error pz ${OpenCV_LIBS})


file(COPY Mazes DESTINATION .)
