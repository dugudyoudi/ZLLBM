# /examples/lbm/immersed_boundary/sphere3d/CMakeLists.txt
set(PROJECT_NAME lbm_ib_sphere3d)
set(EXECUTABLE_NAME "sphere3d")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin/examples/lbm/immersed_boundary/sphere3d)
include_directories(${CMAKE_SOURCE_DIR}/source/lbm_project)
include_directories(${CMAKE_SOURCE_DIR}/source/amr_project)
add_executable(${EXECUTABLE_NAME} "./${EXECUTABLE_NAME}.cpp")
message("${CMAKE_CURRENT_DIR}")
target_link_libraries(${EXECUTABLE_NAME} PRIVATE lbm_project)
