if(BUILD_WITH_PROFILING)
    if(NOT BUILD_WITH_MPI)
        message(FATAL_ERROR "BUILD_WITH_PROFILING requires BUILD_WITH_MPI to be enabled. CMake will quit.")
    endif()

    # Find MPI first (required for MPI profiling)
    find_package(MPI REQUIRED)
    if(NOT MPI_FOUND)
        message(FATAL_ERROR "MPI not found. Please ensure MPI is installed and available.")
    endif()

    # Allow user to specify Score-P root directory via a cache variable
    set(SCOREP_ROOT "" CACHE PATH "Root directory of Score-P installation")

    # Search for Score-P MPI compiler wrappers
    if(SCOREP_ROOT)
        find_program(SCOREP_CXX_COMPILER
            NAMES scorep-mpicxx
            PATHS "${SCOREP_ROOT}/bin"
            NO_DEFAULT_PATH
            DOC "Path to Score-P MPI C++ compiler wrapper"
        )
        find_program(SCOREP_C_COMPILER
            NAMES scorep-mpicc
            PATHS "${SCOREP_ROOT}/bin"
            NO_DEFAULT_PATH
            DOC "Path to Score-P MPI C compiler wrapper"
        )
    else()
        find_program(SCOREP_CXX_COMPILER
            NAMES scorep-mpicxx
            PATHS
                ENV PATH
                /home/$ENV{USER}/scorep/bin
                /opt/scorep/bin
                /usr/local/scorep/bin
                /usr/local/scorep*/bin
                /usr/bin
            DOC "Path to Score-P MPI C++ compiler wrapper"
        )
        find_program(SCOREP_C_COMPILER
            NAMES scorep-mpicc
            PATHS
                ENV PATH
                /home/$ENV{USER}/scorep/bin
                /opt/scorep/bin
                /usr/local/scorep/bin
                /usr/local/scorep*/bin
                /usr/bin
            DOC "Path to Score-P MPI C compiler wrapper"
        )
    endif()
    message(STATUS "Current USER: $ENV{USER}")
    message(STATUS "Current PATH: $ENV{PATH}")
    if(NOT SCOREP_CXX_COMPILER OR NOT SCOREP_C_COMPILER)
        message(STATUS "SCOREP_CXX_COMPILER: ${SCOREP_CXX_COMPILER}")
        message(STATUS "SCOREP_C_COMPILER: ${SCOREP_C_COMPILER}")
        message(FATAL_ERROR "Score-P MPI compiler wrappers (scorep-mpicxx, scorep-mpicc) not found. "
                            "Please install Score-P, add it to PATH, or specify -DSCOREP_ROOT=/path/to/scorep.")
    endif()

    # Check if the wrappers were found
    if(NOT SCOREP_CXX_COMPILER OR NOT SCOREP_C_COMPILER)
        message(FATAL_ERROR "Score-P MPI compiler wrappers (scorep-mpicxx, scorep-mpicc) not found. "
                            "Please install Score-P, add it to PATH, or specify -DSCOREP_ROOT=/path/to/scorep.")
    endif()

    # Override the CMake compilers with Score-P wrappers
    set(CMAKE_CXX_COMPILER "${SCOREP_CXX_COMPILER}" CACHE FILEPATH "C++ compiler overridden by Score-P" FORCE)
    set(CMAKE_C_COMPILER "${SCOREP_C_COMPILER}" CACHE FILEPATH "C compiler overridden by Score-P" FORCE)
    message(STATUS "CMAKE_CXX_COMPILER set to: ${CMAKE_CXX_COMPILER}")
    message(STATUS "CMAKE_C_COMPILER set to: ${CMAKE_C_COMPILER}")

    # Inform user about Score-P configuration
    message(STATUS "Score-P profiling enabled with MPI support")

    # Note: Environment variables like SCOREP_ENABLE_PROFILING are set at runtime, not build time
    # Users need to set these when running the application, e.g., via a script or command line
    message(STATUS "To enable profiling at runtime, set: export SCOREP_ENABLE_PROFILING=1")
    message(STATUS "To enable tracing instead, set: export SCOREP_ENABLE_TRACING=1")
endif()