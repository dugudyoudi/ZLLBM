# /examples/lbm/immersed_boundary/cylinder2d/CMakeLists.txt
set(PROJECT_NAME lbm_ib_cylinder2d)
set(EXECUTABLE_NAME "cylinder2d")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin/examples/lbm/immersed_boundary/cylinder2d)
include_directories(${CMAKE_SOURCE_DIR}/source/lbm_project)
include_directories(${CMAKE_SOURCE_DIR}/source/amr_project)
add_executable(${EXECUTABLE_NAME} "./${EXECUTABLE_NAME}.cpp")
message("${CMAKE_CURRENT_DIR}")
target_link_libraries(${EXECUTABLE_NAME} PRIVATE amr_project lbm_project)
