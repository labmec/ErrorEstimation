add_executable(PerformanceMHM MHM_Performance.cpp)
target_link_libraries(PerformanceMHM Tools ErrorEstimationLib)

add_executable(PerformanceShape shape_Performance.cpp)
target_link_libraries(PerformanceShape NeoPZ::pz)

add_executable(PerformanceContribute Performance\ Contribute.cpp)
target_link_libraries(PerformanceContribute NeoPZ::pz)

add_executable(PerformanceComputeRequiredData Performance\ ComputeRequiredData.cpp)
target_link_libraries(PerformanceComputeRequiredData NeoPZ::pz)
