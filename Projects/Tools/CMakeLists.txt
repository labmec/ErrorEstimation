# This file creates a library with tools files
add_library(Tools
    Tools.cpp
    Tools.h
    ToolsMHM.cpp
    ToolsMHM.h
    DataStructure.h
)

target_link_libraries(Tools ErrorEstimationLib NeoPZ::pz)
target_include_directories(Tools PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
