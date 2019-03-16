# Project to test Mark Ainsworth purpose 
file(GLOB header *.h)
file(GLOB source *.cpp)

set(ErrorEstimateSRC ${header} ${source})

add_executable(PropostaMark ${ErrorEstimateSRC})

target_link_libraries(PropostaMark pz)

configure_file(LMesh.msh LMesh.msh COPYONLY)
