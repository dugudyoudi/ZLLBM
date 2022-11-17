cmake_minimum_required(VERSION 3.1)

# project information
set(PROJECT_NAME rootProject)
set(PROJECT_VERSION 0.1.0)

include(${CMAKE_SOURCE_DIR}/cmake/build_functions.cmake)

# settings
option(BUILD_WITH_MPI                          "Enable MPI"                                                OFF)
option(BUILD_TEST_CASES                        "Build cases for test"                                      ON)

option(DEBUG_UNIT_TEST                         "Buil unit tests"                                           ON)
option(DEBUG_CHECK_GRID                        "Enable if clauses for checking grid"                       ON)
option(DEBUG_BUILD_2D_ONLY                     "Disable 3D related functions"                              OFF)
option(DEBUG_BUILD_3D_ONLY                     "Disable 2D related functions"                              OFF)
option(DEBUG_TREAT_WARNING_AS_ERROR            "Treat warnings as errors"                                  ON)

set(CMAKE_BUILD_TYPE "Debug")
set(CXX_STANDARD "c++20")

# config project
project(${PROJECT_NAME}
    VERSION ${PROJECT_VERSION}
)

enable_language(CXX)
if(BUILD_WITH_MPI)
    set(CONFIG_MPI "#define ENABLE_MPI")
endif(BUILD_WITH_MPI)
if(DEBUG_CHECK_GRID)
    set(CONFIG_DEBUG_CHECK_GRID "#define DEBUG_CHECK_GRID")
endif(DEBUG_CHECK_GRID)
if(DEBUG_BUILD_2D_ONLY )
    # define command disabling 3D related functions
    set(CONFIG_DEBUG_BUILD_2D_ONLY "#define DEBUG_DISABLE_3D_FUNCTIONS")
    if (DEBUG_BUILD_3D_ONLY)
        message(FATAL_ERROR "Can't disable 2D and 3D functions simultaneously. CMake will quit.")
    endif(DEBUG_BUILD_3D_ONLY)
endif(DEBUG_BUILD_2D_ONLY)
if(DEBUG_BUILD_3D_ONLY)
    # define command disabling 2D related functions
    set(CONFIG_DEBUG_BUILD_3D_ONLY "#define DEBUG_DISABLE_2D_FUNCTIONS")
endif(DEBUG_BUILD_3D_ONLY)
if(DEBUG_UNIT_TEST)
    # define command for unit test to declare friend class
    set(CONFIG_DEBUG_UNIT_TEST "#define DEBUG_UNIT_TEST")
endif(DEBUG_UNIT_TEST)
configure_file (
"${CMAKE_SOURCE_DIR}/config.h.in"
"${CMAKE_SOURCE_DIR}/source/config.h"
)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()
message("Build Type: ${CMAKE_BUILD_TYPE}")

# System
if (CMAKE_SYSTEM_NAME MATCHES "Linux")
    set(CURRENT_SYSTEM "Linux")
elseif (CMAKE_SYSTEM_NAME MATCHES "Windows")
    set(CURRENT_SYSTEM "Windows")
else()
    set(CURRENT_SYSTEM "Other system")
endif()

if(CMAKE_CL_64)
    set(CURRENT_PLATFORM "x64")
else(CMAKE_CL_64)
    set(CURRENT_PLATFORM "x86")
endif(CMAKE_CL_64)
message(STATUS "Current Platform is ${CURRENT_SYSTEM} ${CURRENT_PLATFORM}")

# Compiler
if(${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
    message("Compiler: MSVC")
    set(CMAKE_CXX_FLAGS "/std:${CXX_STANDARD}")
    set(CMAKE_CXX_FLAGS_RELEASE "/Ox /Ot /GS- /GL")
    set(CMAKE_CXX_FLAGS_DEBUG "/W3 /Zi /EHsc")
    if (DEBUG_TREAT_WARNING_AS_ERROR)
      AppendFlag(CMAKE_CXX_FLAGS_DEBUG "/WX")
    endif()
elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL GNU)
    message("Compiler: GNU")
    set(CMAKE_CXX_FLAGS "-std=${CXX_STANDARD} -Wall -Wextra -Wwrite-strings -Wno-parentheses -Wpedantic -Warray-bounds  -Weffc++")
    set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -fexec-charset=GBK")
    set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
    if (DEBUG_TREAT_WARNING_AS_ERROR)
      AppendFlag(CMAKE_CXX_FLAGS_DEBUG "-pedantic-errors -Werror")
      message(CMAKE_CXX_FLAGS_DEBUG)
    endif()
elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL Clang)
    message("Compiler: Clang")
    set(CMAKE_CXX_FLAGS "-std=${CXX_STANDARD} -Wall -Wextra -Wwrite-strings -Wno-parentheses -Wpedantic -Warray-bounds  -Weffc++")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3")
    set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
    if (DEBUG_TREAT_WARNING_AS_ERROR)
      AppendFlag(CMAKE_CXX_FLAGS_DEBUG "-pedantic-errors -Werror")
    endif()
else()
    message(FATAL_ERROR "CXX compiler not recognized. CMake will quit.")
endif()

# MSMPI is chosen priorly 
if(BUILD_WITH_MPI)
    find_package(MPI REQUIRED)
    if(MPI_CXX_FOUND)
        if(EXISTS "$ENV{MSMPI_INC}") # visual studio default mpi is impi, can't initialize correctly.
           include_directories("$ENV{MSMPI_INC}"  "$ENV{MSMPI_INC}/${CURRENT_PLATFORM}")
           message("MSMPI directories: $ENV{MSMPI_INC}")
        else(DEFINED $ENV{MSMPI_INC})
           include_directories(${MPI_CXX_INCLUDE_PATH})
           message("${MPI_CXX_COMPILER} directories: ${MPI_CXX_INCLUDE_DIRS}")  
        endif(EXISTS "$ENV{MSMPI_INC}")

        set(CMAKE_CXX_FLAGS"${CMAKE_CXX_FLAGS} ${MPI_CXX_FLAGS}")
    endif(MPI_CXX_FOUND)
endif(BUILD_WITH_MPI)

add_subdirectory("./source/amr_project")

# test cases
if(BUILD_TEST_CASES)
    add_subdirectory("test_cases/test_grid")
endif(BUILD_TEST_CASES)

# unit test
if(DEBUG_UNIT_TEST)
    # compile gtest
    add_subdirectory("./googletest-main")
    # compile test cases
    add_subdirectory("./unit_tests/test_amr")
endif(DEBUG_UNIT_TEST)



