# common source files used by targets of project ErrorEstimation

# enabling to use LOG4CXX library
option(USING_LOG4CXX "Whether the LOG4CXX library will be linked in" ON)
if(USING_LOG4CXX)
    find_package(log4cxx)
endif(USING_LOG4CXX)

# enabling PZ library
IF (WIN32)
    # Prevents timespec redefinition problem with Visual Studio 2015
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHAVE_STRUCT_TIMESPEC")
    find_package( PZ PATHS "C:/Arquivos de Programas/PZ" REQUIRED NO_DEFAULT_PATH)
ELSE ()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-narrowing")#this flag is for preventing errors when casting from double to float
    find_package( PZ REQUIRED )
ENDIF()

IF(APPLE)
    # Mac OS X specific code
    add_definitions(-DMACOSX)
    set(MACOSX MACOSX)
    FIND_LIBRARY(ACCELERATE_LIB Accelerate)
    link_libraries(${ACCELERATE_LIB})
    #   SET(OperatingSystem "Mac OS X")
ENDIF(APPLE)

if(WIN32)
    #add a define saying it is a VC compiler
    set(VC "#Will add a define VC on pz_config.h")
    add_definitions(-DVC)
    #define _USE_MATH_DEFINES for example: M_PI constant definitions.
    add_definitions(-D_USE_MATH_DEFINES)
    #disabling VC warnings
    add_definitions(-D_SCL_SECURE_NO_WARNINGS)
    add_definitions(-D_CRT_SECURE_NO_WARNINGS)
    #disabling All warnings
    #set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /W0")
    #set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /W0")
    #disabling STL warning for std::_Vector_alloc when LOG4CXX is enabled
    if (USING_LOG4CXX)
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /wd4251")
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /wd4251")
    endif()
    #define use of pthread static lib.
    #enabling /bigobj
    add_definitions("/bigobj")

    FOREACH(FLAG_TYPE EXE MODULE SHARED)
        STRING (REPLACE "INCREMENTAL:YES" "INCREMENTAL" FLAG_TMP "${CMAKE_${FLAG_TYPE}_LINKER_FLAGS_DEBUG}")
        STRING (REPLACE "INCREMENTAL:NO" "INCREMENTAL" FLAG_TMP ${FLAG_TMP})
        STRING (REPLACE "INCREMENTAL" "INCREMENTAL:NO" FLAG_TMP ${FLAG_TMP})
        STRING (REPLACE "/EDITANDCONTINUE" "" FLAG_TMP ${FLAG_TMP})
        SET(CMAKE_${FLAG_TYPE}_LINKER_FLAGS_DEBUG "${FLAG_TMP}" CACHE STRING "Overriding default debug ${FLAG_TYPE} linker flags." FORCE)
        MARK_AS_ADVANCED (CMAKE_${FLAG_TYPE}_LINKER_FLAGS_DEBUG)
    ENDFOREACH ()

endif()


file(GLOB all_sources *.cpp *.h)
set(common_SRC ${all_sources} PARENT_SCOPE)

add_library (Common ${all_sources})
target_include_directories (Common PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${PZ_INCLUDE_DIRS})
target_link_libraries(Common PUBLIC pz)
