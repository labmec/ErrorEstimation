

set(MHMMAZE mainMaze.cpp ConfigCasesMaze.h meshgen.h meshgen.cpp)

add_executable(MHM_Maze mainMaze.cpp ${MHMMAZE})
target_link_libraries(MHM_Maze pz ${OpenCV_LIBS})

