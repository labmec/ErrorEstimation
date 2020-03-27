add_executable(UNISIMHDiv mainHdiv.cpp ${common_SRC})

find_package(libInterpolate REQUIRED)
target_link_libraries(UNISIMHDiv
        PUBLIC libInterpolate::Interpolate pz
        )

if (ENABLE_ALL_WARNINGS)
    target_link_libraries(UNISIMHDiv PRIVATE project_warnings)
endif ()

configure_file(InputData/UNISIMFlatMesh.msh InputData/UNISIMFlatMesh.msh COPYONLY)
configure_file(InputData/UNISIMPointCloud.txt InputData/UNISIMPointCloud.txt COPYONLY)
