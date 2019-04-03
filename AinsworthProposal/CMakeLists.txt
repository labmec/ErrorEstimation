# This target implements M. Ainsworth proposal using new NeoPZ functionalities.

file(GLOB header *.h)
file(GLOB source *.cpp)

set(ErrorEstimateSRC ${header} ${source})

add_executable(AinsworthProposal ${ErrorEstimateSRC})

target_link_libraries(AinsworthProposal pz)

configure_file(LMesh.msh LMesh.msh COPYONLY)
