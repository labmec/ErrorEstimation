find_package(OpenCV REQUIRED)

add_executable(MHM_Maze main_Maze.cpp Maze_common.cpp Maze_common.h ConfigCasesMaze.h)
target_link_libraries(MHM_Maze Tools ErrorEstimationLib ${OpenCV_LIBS})

add_executable(MHM_Steklov main_Maze_Steklov.cpp Maze_common.cpp Maze_common.h ConfigCasesMaze.h)
target_link_libraries(MHM_Steklov Tools ErrorEstimationLib ${OpenCV_LIBS})

add_executable(MHM_Steklov2 main_Maze_Steklov_constant.cpp Maze_common.cpp Maze_common.h ConfigCasesMaze.h)
target_link_libraries(MHM_Steklov2 Tools ErrorEstimationLib ${OpenCV_LIBS})

find_package(libInterpolate REQUIRED)
add_executable(spe10_Steklov main_spe10.cpp Maze_common.cpp Maze_common.h ConfigCasesMaze.h)
target_link_libraries(spe10_Steklov Tools ErrorEstimationLib ${OpenCV_LIBS} libInterpolate::Interpolate)

file(COPY Mazes DESTINATION .)
