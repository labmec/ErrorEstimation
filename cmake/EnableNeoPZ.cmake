function(enable_pz)

    find_package(PZ REQUIRED)

    include_directories(${PZ_INCLUDE_DIRS})

endfunction()
