add_executable(UNISIMHDiv mainHdiv.cpp ${common_SRC} )

target_link_libraries(UNISIMHDiv pz)

find_package(libInterpolate REQUIRED)
target_link_libraries(UNISIMHDiv libInterpolate::Interpolate)

configure_file(InputData/UNISIMFlatMesh.msh InputData/UNISIMFlatMesh.msh COPYONLY)
configure_file(InputData/UNISIMPointCloud.txt InputData/UNISIMPointCloud.txt COPYONLY)
