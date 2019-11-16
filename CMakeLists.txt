# @file SBFem/CMakeList.txt  -- First file to read for CMake

cmake_minimum_required (VERSION 3.10.0)

project (ErrorEstimation)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# Creates library with common files
add_subdirectory(Common)

# Projects
add_subdirectory(ErrorEstimationH1)
add_subdirectory(ErrorEstimationHDiv)
add_subdirectory(HybridH1)
add_subdirectory(MHMEstimation)

# OpenCV projects
option(BUILD_OPENCV_PROJECTS "Whether to build projects that require OpenCV" OFF)
if(BUILD_OPENCV_PROJECTS)
	add_subdirectory(MHM_Maze) # must be added here to inherit non-empty ${common_SRC}
endif()

# Python bindings
option(BUILD_PYTHON_BINDINGS "Whether to build python bindings" ON)
if(BUILD_PYTHON_BINDINGS)
	add_subdirectory(python-bindings)
endif()
