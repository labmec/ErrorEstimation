function(enable_pz)

IF (WIN32)
    # Prevents timespec redefinition problem with Visual Studio 2015
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHAVE_STRUCT_TIMESPEC")
    find_package(PZ PATHS "C:/Arquivos de Programas/PZ" REQUIRED NO_DEFAULT_PATH)
ELSE ()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-narrowing")#this flag is for preventing errors when casting from double to float
    find_package(PZ REQUIRED)
ENDIF ()

IF(APPLE)
	# Mac OS X specific code
	add_definitions(-DMACOSX)
	set(MACOSX MACOSX)
	find_library(ACCELERATE_LIB Accelerate)
	link_libraries(${ACCELERATE_LIB})
	#   SET(OperatingSystem "Mac OS X")
ENDIF(APPLE)

include_directories(${PZ_INCLUDE_DIRS})

endfunction()
