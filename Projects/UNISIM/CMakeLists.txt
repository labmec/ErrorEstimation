# Hybrid Hdiv UNISIM simulation target
add_executable(UNISIM_HDiv mainHdiv.cpp ToolsUNISIM.cpp)

find_package(libInterpolate REQUIRED)
target_link_libraries(UNISIM_HDiv PUBLIC
    libInterpolate::Interpolate
    Tools
    ErrorEstimationLib
)

# MHM-Hdiv UNISIM simulation target
add_executable(UNISIM_MHM mainMHM.cpp ToolsUNISIM.cpp)

find_package(libInterpolate REQUIRED)
target_link_libraries(UNISIM_MHM PUBLIC
    libInterpolate::Interpolate
    Tools
    ErrorEstimationLib
    )

configure_file(InputData/UNISIMFlatMesh.msh InputData/UNISIMFlatMesh.msh COPYONLY)
configure_file(InputData/UNISIMPointCloud.txt InputData/UNISIMPointCloud.txt COPYONLY)
