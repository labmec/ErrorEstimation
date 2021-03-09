cmake_minimum_required (VERSION 3.13)

project (ErrorEstimation)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

include(cmake/StandardProjectSettings.cmake)

# Enables NeoPZ library
include(cmake/EnableNeoPZ.cmake)
enable_pz()

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

option(USING_LOG4CXX "Whether the LOG4CXX library will be linked in" ON)
if(USING_LOG4CXX)
    #this is the definition used in pz lib
    find_package(Log4cxx REQUIRED)
    target_link_libraries(ErrorEstimationLib PRIVATE ${Log4cxx_LIBRARY})
    target_include_directories(ErrorEstimationLib PRIVATE ${Log4cxx_INCLUDE_DIR})
endif(USING_LOG4CXX)
