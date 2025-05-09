# /lbm_projectrce
cmake_minimum_required(VERSION 3.15...3.5)

# project information
set(LBM_PROJECT_NAME lbm_project)
set(PROJECT_VERSION 0.1.0)

# config project
project(${LBM_PROJECT_NAME}
    VERSION ${PROJECT_VERSION}
)

# directories of libs and sources
include_directories("${CMAKE_CURRENT_SOURCE_DIR}")
FILE(GLOB SOURCE_LBM_GENERAL_SRC "./*.cpp")
FILE(GLOB SOURCE_LBM_GENERAL_HEAD "./*.h")
FILE(GLOB SOURCE_LBM_BOUNDARY_SRC "./boundary_conditions/*.cpp")
FILE(GLOB SOURCE_LBM_BOUNDARY_HEAD "./boundary_conditions/*.h")
FILE(GLOB SOURCE_LBM_IMMERSED_BOUNDARY_SRC "./immersed_boundary/*.cpp")
FILE(GLOB SOURCE_LBM_IMMERSED_BOUNDARY_HEAD "./immersed_boundary/*.h")
FILE(GLOB SOURCE_LBM_COLLISION_MODEL_SRC "./collision_model/*.cpp")
FILE(GLOB SOURCE_LBM_COLLISION_MODEL_HEAD "./collision_model/*.h")
FILE(GLOB SOURCE_LBM_LES_MODEL_SRC "./les_models/*.cpp")
FILE(GLOB SOURCE_LBM_LES_MODEL_HEAD "./les_models/*.h")
FILE(GLOB SOURCE_LBM_D2Q9_SRC "./d2q9_model/*.cpp")
FILE(GLOB SOURCE_LBM_D2Q9_HEAD "./d2q9_model/*.h")
FILE(GLOB SOURCE_LBM_D3Q19_SRC "./d3q19_model/*.cpp")
FILE(GLOB SOURCE_LBM_D3Q19_HEAD "./d3q19_model/*.h")

add_library(${LBM_PROJECT_NAME} STATIC 
            ${SOURCE_LBM_GENERAL_SRC} ${SOURCE_LBM_GENERAL_HEAD}
            ${SOURCE_LBM_BOUNDARY_SRC} ${SOURCE_LBM_BOUNDARY_HEAD}
            ${SOURCE_LBM_IMMERSED_BOUNDARY_SRC} ${SOURCE_LBM_IMMERSED_BOUNDARY_HEAD}
            ${SOURCE_LBM_COLLISION_MODEL_SRC} ${SOURCE_LBM_COLLISION_MODEL_HEAD}
            ${SOURCE_LBM_LES_MODEL_SRC} ${SOURCE_LBM_LES_MODEL_HEAD}
            ${SOURCE_LBM_D2Q9_SRC} ${SOURCE_LBM_D2Q9_HEAD}
            ${SOURCE_LBM_D3Q19_SRC} ${SOURCE_LBM_D3Q19_HEAD})

target_link_libraries(${LBM_PROJECT_NAME} amr_project)
if(BUILD_WITH_MPI)
    target_link_libraries(${LBM_PROJECT_NAME} ${BUILD_LINK_LIBRARY})
endif(BUILD_WITH_MPI)
if(WIN32)
    target_link_libraries(${LBM_PROJECT_NAME} ws2_32)
endif(WIN32)

# setup solution filter
source_group("HeadFiles" FILES ${SOURCE_LBM_HEAD})
source_group("SourceFiles" FILES ${SOURCE_LBM_SRC})
source_group("HeadFiles\\boundary_conditions" FILES ${SOURCE_LBM_BOUNDARY_HEAD})
source_group("SourceFiles\\boundary_conditions"  FILES ${SOURCE_LBM_BOUNDARY_SRC})
source_group("HeadFiles\\immersed_boundary" FILES ${SOURCE_LBM_IMMERSED_BOUNDARY_HEAD})
source_group("SourceFiles\\immersed_boundary"  FILES ${SOURCE_LBM_IMMERSED_BOUNDARY_SRC})
source_group("HeadFiles\\collision_model" FILES ${SOURCE_LBM_COLLISION_MODEL_HEAD})
source_group("SourceFiles\\collision_model"  FILES ${SOURCE_LBM_COLLISION_MODEL_SRC})
source_group("HeadFiles\\les_models" FILES ${SOURCE_LBM_LES_MODEL_HEAD})
source_group("SourceFiles\\les_models"  FILES ${SOURCE_LBM_LES_MODEL_SRC})
source_group("HeadFiles\\d2q9_model"  FILES ${SOURCE_LBM_D2Q9_HEAD})
source_group("SourceFiles\\d2q9_model" FILES ${SOURCE_LBM_D2Q9_SRC})
source_group("HeadFiles\\d3q19_model"  FILES ${SOURCE_LBM_D3Q19_HEAD})
source_group("SourceFiles\\d3q19_model" FILES ${SOURCE_LBM_D3Q19_SRC})

# include amr project
include_directories(${CMAKE_SOURCE_DIR}/source/amr_project)
# include unit test
if(DEBUG_UNIT_TEST)
    include_directories(${CMAKE_SOURCE_DIR}/googletest-main/googletest/include)
    include_directories(${CMAKE_SOURCE_DIR}/googletest-main/googlemock/include)
endif(DEBUG_UNIT_TEST)