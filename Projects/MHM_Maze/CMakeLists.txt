find_package(OpenCV REQUIRED)

add_executable(MHM_Maze mainMaze.cpp ConfigCasesMaze.h)
target_link_libraries(MHM_Maze Tools ErrorEstimationLib ${OpenCV_LIBS})

add_executable(MHM_Error mainMHMError.cpp ConfigCasesMaze.h)
target_link_libraries(MHM_Error Tools ErrorEstimationLib ${OpenCV_LIBS})

file(COPY Mazes DESTINATION .)
