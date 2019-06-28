

set(MHMMAZE mainMaze.cpp ConfigCasesMaze.h meshgen.h meshgen.cpp
)

add_executable(MHM_Maze mainMaze.cpp ${MHMMAZE})
target_link_libraries(MHM_Maze pz ${OpenCV_LIBS})

set(MHMError mainMHMError.cpp ConfigCasesMaze.h meshgen.h meshgen.cpp ${common_SRC})

add_executable(MHM_Error ${MHMError})
target_link_libraries(MHM_Error pz ${OpenCV_LIBS})


file(COPY Mazes DESTINATION .)
