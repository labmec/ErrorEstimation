find_package(OpenCV REQUIRED)

add_executable(MHM_Maze main_Maze.cpp ConfigCasesMaze.h)
target_link_libraries(MHM_Maze Tools ErrorEstimationLib ${OpenCV_LIBS})

file(COPY Mazes DESTINATION .)
