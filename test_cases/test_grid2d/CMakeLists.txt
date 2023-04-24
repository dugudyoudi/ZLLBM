# /test_cases/test_grid/CMakeLists.txt
set(PROJECT_NAME test_cases)
set(PROJECT_VERSION "0.1.0")
set(EXECUTABLE_NAME "test_grid2d")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin/grid_test)
include_directories(${CMAKE_SOURCE_DIR}/source/amr_project)
add_executable(${EXECUTABLE_NAME} "./${EXECUTABLE_NAME}.cpp")
target_link_libraries(${EXECUTABLE_NAME} amr_project)
