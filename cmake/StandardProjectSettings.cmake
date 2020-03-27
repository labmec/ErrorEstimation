# Enables ccache to increase compilation speed
option(ENABLE_CCACHE "Enable ccache" OFF)
if (ENABLE_CCACHE)
    find_program(CCACHE ccache)
    if (CCACHE)
        message("using ccache")
        set(CMAKE_CXX_COMPILER_LAUNCHER ${CCACHE})
    else ()
        message("ccache not found cannot use")
    endif ()
endif ()

# Generate compile_commands.json to make it easier to work with clang based
# tools
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(ENABLE_IPO "Enable Interprocedural Optimization" OFF)
if (ENABLE_IPO)
    include(CheckIPOSupported)
    check_ipo_supported(RESULT result OUTPUT output)
    if (result)
        set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
    else ()
        message(SEND_ERROR "IPO is not supported: ${output}")
    endif ()
endif ()
