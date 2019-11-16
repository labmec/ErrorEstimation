#ckwg +4
# Copyright 2010 by Kitware, Inc. All Rights Reserved. Please refer to
# KITWARE_LICENSE.TXT for licensing information, or contact General Counsel,
# Kitware, Inc., 28 Corporate Drive, Clifton Park, NY 12065.

# Locate the system installed log4cxx
# The following variables will be set:
#
# log4cxx_FOUND       - Set to true if log4cxx can be found
# log4cxx_INCLUDE_DIR - The path to the log4cxx header files
# log4cxx_LIBRARY     - The full path to the log4cxx library

if( log4cxx_DIR )
    find_package( log4cxx NO_MODULE )
elseif( NOT log4cxx_FOUND )
    message(STATUS "Searching for log4cxx/logger.h")
    find_path( log4cxx_INCLUDE_DIR log4cxx/logger.h )

    message(STATUS "Searching for liblog4cxx")
    find_library( log4cxx_LIBRARY log4cxx )

    include( FindPackageHandleStandardArgs )
    FIND_PACKAGE_HANDLE_STANDARD_ARGS( log4cxx log4cxx_INCLUDE_DIR log4cxx_LIBRARY )
    if( log4cxx_FOUND )
        set( log4cxx_FOUND TRUE )
        mark_as_advanced(log4cxx_INCLUDE_DIR log4cxx_LIBRARY)
    endif()
endif()