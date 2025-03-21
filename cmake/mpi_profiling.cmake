# Only proceed if BUILD_WITH_TAU is enabled
if(BUILD_WITH_TAU)
    if(NOT BUILD_WITH_MPI)
        message(FATAL_ERROR "BUILD_WITH_TAU requires BUILD_WITH_MPI to be enabled. CMake will quit.")
    endif()

    # Allow user to specify TAU root directory via cache variable
    set(TAU_ROOT "" CACHE PATH "Root directory of TAU installation")
    if(TAU_ROOT)
        find_program(TAU_CXX_COMPILER
            NAMES tau_cxx.sh
            PATHS "${TAU_ROOT}/x86_64/bin"
            NO_DEFAULT_PATH
            DOC "Path to TAU C++ compiler wrapper"
        )
    else()
        find_program(TAU_CXX_COMPILER
            NAMES tau_cxx.sh
            PATHS
                ENV PATH
                /home/${USER}/tau*/x86_64/bin
                /opt/tau*/x86_64/bin
                /usr/local/tau*/x86_64/bin
            DOC "Path to TAU C++ compiler wrapper"
        )
    endif()
    if(NOT TAU_CXX_COMPILER)
        message(FATAL_ERROR "TAU compiler (tau_cxx.sh) not found. Please install TAU, add it to PATH, or specify -DTAU_ROOT=/path/to/tau.")
    endif()

    # Override the C++ compiler with TAU's wrapper
    set(CMAKE_CXX_COMPILER "${TAU_CXX_COMPILER}" CACHE FILEPATH "C++ compiler overridden by TAU" FORCE)

    # Set TAU environment variables
    get_filename_component(TAU_BIN_DIR "${CMAKE_CXX_COMPILER}" DIRECTORY)
    get_filename_component(TAU_BASE_DIR "${TAU_BIN_DIR}" DIRECTORY)
    set(TAU_MAKEFILE_PATH "${TAU_BASE_DIR}/lib/Makefile.tau-mpi")
    if(EXISTS "${TAU_MAKEFILE_PATH}")
        set(ENV{TAU_MAKEFILE} "${TAU_MAKEFILE_PATH}")
    else()
        message(FATAL_ERROR "TAU Makefile not found at ${TAU_MAKEFILE_PATH}. Check TAU installation.")
    endif()
    message(STATUS "TAU profiling enabled with compiler: ${CMAKE_CXX_COMPILER}")
endif()