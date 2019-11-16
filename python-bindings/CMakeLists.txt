# pybind11
find_package(pybind11 REQUIRED)

# Project related sources
pybind11_add_module(errorestimation src/bindings.cpp)
target_link_libraries(errorestimation PUBLIC Common)
##target_include_directories(errorestimationbindings PUBLIC ${INCLUDE_DIR})
