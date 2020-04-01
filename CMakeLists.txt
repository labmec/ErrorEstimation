cmake_minimum_required (VERSION 3.13)

project (ErrorEstimation)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(cmake/StandardProjectSettings.cmake)

## Sanitizer options if supported by compiler
#include(cmake/Sanitizers.cmake)
#add_library(project_options INTERFACE)
#enable_sanitizers(project_options)

# Allow for static analysis options
include(cmake/StaticAnalyzers.cmake)

# Enables NeoPZ library
include(cmake/EnablePZ.cmake)
enable_pz()

# Enables log4cxx library
include(cmake/EnableLog4CXX.cmake)
enable_log4cxx()

# This option enables a lot of warnings and treat them as errors, to ensure
# good programming practices are used. Since its behaviour is extreme, it
# should be turned off by default.
option(ENABLE_ALL_WARNINGS "Enable a lot of compiler warnings and treat them as errors" OFF)
if (ENABLE_ALL_WARNINGS)
    add_library(project_warnings INTERFACE)
    include(cmake/CompilerWarnings.cmake)
    set_project_warnings(project_warnings)
    set(CMAKE_CXX_FLAGS "")
endif ()

add_subdirectory(Common)
include_directories(Common)
add_subdirectory(ErrorEstimationH1)
add_subdirectory(ErrorEstimationHDiv)
add_subdirectory(HybridH1)
add_subdirectory(MHMEstimation)

# This project is marked as optional since it requires OpenCV external library
option(BUILD_MHM_MAZE_PROJECT "Whether to build MHM Maze project" OFF)
if(BUILD_MHM_MAZE_PROJECT)
    add_subdirectory(MHM_Maze)
endif()

# This project is marked as optional since it requires libInterpolate external library
option(BUILD_MHM_MAZE_PROJECT "Whether to build MHM Maze project" OFF)
option(BUILD_UNISIM_PROJECT "Whether to build estimation on UNISIM mesh" OFF)
if(BUILD_UNISIM_PROJECT)
    add_subdirectory(UNISIM)
endif()
