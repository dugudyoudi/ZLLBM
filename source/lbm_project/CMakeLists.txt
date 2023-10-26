# /source

cmake_minimum_required(VERSION 3.5)

# project information
set(LBM_PROJECT_NAME lbm_project)
set(PROJECT_VERSION 0.1.0)

# config project
project(${LBM_PROJECT_NAME}
    VERSION ${PROJECT_VERSION}
)

# directories of libs and sources
include_directories("${CMAKE_CURRENT_SOURCE_DIR}")
FILE(GLOB SOURCE_LBM_SRC "./*.cpp")
FILE(GLOB SOURCE_LBM_HEAD "./*.h")

add_library(${LBM_PROJECT_NAME} STATIC 
            ${SOURCE_LBM_SRC} ${SOURCE_LBM_HEAD})

if(BUILD_WITH_MPI)
    target_link_libraries(${LBM_PROJECT_NAME}
    ${BUILD_LINK_LIBRARY})
endif(BUILD_WITH_MPI)
if(WIN32)
    target_link_libraries(${LBM_PROJECT_NAME} ws2_32)
endif(WIN32)

# setup solution filter
source_group("HeadFiles" FILES ${SOURCE_LBM_HEAD})
source_group("SourceFiles" FILES ${SOURCE_LBM_SRC})

# include amr project
include_directories(${CMAKE_SOURCE_DIR}/source/amr_project)
# include unit test
if(DEBUG_UNIT_TEST)
    include_directories(${CMAKE_SOURCE_DIR}/googletest-main/googletest/include)
    include_directories(${CMAKE_SOURCE_DIR}/googletest-main/googlemock/include)
endif(DEBUG_UNIT_TEST)