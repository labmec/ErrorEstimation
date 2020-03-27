set(MHMMAZE mainMaze.cpp ConfigCasesMaze.h )

add_executable(MHM_Maze mainMaze.cpp ${common_SRC} ${MHMMAZE})

find_package(OpenCV REQUIRED)
target_link_libraries(MHM_Maze pz ${OpenCV_LIBS})

set(MHMError mainMHMError.cpp ConfigCasesMaze.h ${common_SRC})

add_executable(MHM_Error ${common_SRC} ${MHMError})
target_link_libraries(MHM_Error pz ${OpenCV_LIBS})

file(COPY Mazes DESTINATION .)
