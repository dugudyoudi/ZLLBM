# /source

cmake_minimum_required(VERSION 3.1)

# project information
set(AMR_PROJECT_NAME amr_project)
set(PROJECT_VERSION 0.1.0)

# config project
project(${AMR_PROJECT_NAME}
    VERSION ${PROJECT_VERSION}
)

# directories of libs and sources
include_directories("${CMAKE_CURRENT_SOURCE_DIR}")
FILE(GLOB SOURCE_AMR_SRC "./*.cpp")
FILE(GLOB SOURCE_AMR_HEAD "./*.h")
FILE(GLOB SOURCE_AMR_GRID_SRC "./grid/*.cpp")
FILE(GLOB SOURCE_AMR_GRID_HEAD "./grid/*.h" "./grid/*.hpp")
FILE(GLOB SOURCE_AMR_IO_SRC "./io/*.cpp")
FILE(GLOB SOURCE_AMR_IO_HEAD "./io/*.h")
FILE(GLOB SOURCE_AMR_GEOMETRY_SRC "./criterion/*.cpp")
FILE(GLOB SOURCE_AMR_GEOMETRY_HEAD "./criterion/*.h")
FILE(GLOB SOURCE_AMR_MPI_SRC "./mpi/*.cpp")
FILE(GLOB SOURCE_AMR_MPI_HEAD "./mpi/*.h")

add_library(${AMR_PROJECT_NAME} STATIC 
            ${SOURCE_AMR_SRC} ${SOURCE_AMR_HEAD} 
            ${SOURCE_AMR_GRID_SRC} ${SOURCE_AMR_GRID_HEAD}
            ${SOURCE_AMR_IO_SRC} ${SOURCE_AMR_IO_HEAD}
            ${SOURCE_AMR_GEOMETRY_SRC} ${SOURCE_AMR_GEOMETRY_HEAD}
            ${SOURCE_AMR_MPI_SRC} ${SOURCE_AMR_MPI_HEAD})

if(BUILD_WITH_MPI)
    target_link_libraries(${AMR_PROJECT_NAME}
    ${BUILD_LINK_LIBRARY})
endif(BUILD_WITH_MPI)

# setup solution filter
source_group("HeadFiles" FILES ${SOURCE_AMR_HEAD})
source_group("SourceFiles" FILES ${SOURCE_AMR_SRC})
source_group("HeadFiles\\grid" FILES ${SOURCE_AMR_GRID_HEAD})
source_group("SourceFiles\\grid"  FILES ${SOURCE_AMR_GRID_SRC})
source_group("HeadFiles\\io"  FILES ${SOURCE_AMR_IO_HEAD})
source_group("SourceFiles\\io" FILES ${SOURCE_AMR_IO_SRC})
source_group("HeadFiles\\criterion" FILES ${SOURCE_AMR_GEOMETRY_HEAD})
source_group("SourceFiles\\criterion" FILES ${SOURCE_AMR_GEOMETRY_SRC})
source_group("HeadFiles\\mpi" FILES ${SOURCE_AMR_MPI_HEAD})
source_group("SourceFiles\\mpi" FILES ${SOURCE_AMR_MPI_SRC})

# unit test
if(DEBUG_UNIT_TEST)
    include_directories(${CMAKE_SOURCE_DIR}/googletest-main/googletest/include)
    include_directories(${CMAKE_SOURCE_DIR}/googletest-main/googlemock/include)
endif(DEBUG_UNIT_TEST)