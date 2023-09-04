cmake_minimum_required (VERSION 3.13)

project (ErrorEstimation)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

include(cmake/StandardProjectSettings.cmake)

# Enables NeoPZ library
find_package(NeoPZ REQUIRED)

# Mac OS X specific code
if (APPLE)
  add_definitions(-DMACOSX)
endif()


#TODO: Define this using modern CMake standards. @Gustavo, do you know how to do this?
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    option(ERRORESTIMATION_DEBUG "Debug verifications" OFF)
    if(ERRORESTIMATION_DEBUG)
        add_definitions(-DERRORESTIMATION_DEBUG)
    endif()
endif()

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

# Creates ErrorEstimationLib with common error estimation file
add_subdirectory(ErrorEstimation)

# Created executable targets
add_subdirectory(Projects)
