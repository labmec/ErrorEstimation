find_package(OpenCV REQUIRED)

add_executable(MHM_Maze main_Maze.cpp Maze_common.cpp Maze_common.h ConfigCasesMaze.h)
target_link_libraries(MHM_Maze Tools ErrorEstimationLib ${OpenCV_LIBS})

add_executable(MHM_Steklov main_Maze_Steklov.cpp Maze_common.cpp Maze_common.h ConfigCasesMaze.h)
target_link_libraries(MHM_Steklov Tools ErrorEstimationLib ${OpenCV_LIBS})


file(COPY Mazes DESTINATION .)
