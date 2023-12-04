# /run_cases/test_tmp/CMakeLists.txt
set(PROJECT_NAME test_cases)
set(PROJECT_VERSION "0.1.0")
set(EXECUTABLE_NAME "test_tmp")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin/run_cases/test_tmp)
include_directories(${CMAKE_SOURCE_DIR}/source/lbm_project)
include_directories(${CMAKE_SOURCE_DIR}/source/amr_project)
add_executable(${EXECUTABLE_NAME} "./${EXECUTABLE_NAME}.cpp")
target_link_libraries(${EXECUTABLE_NAME} amr_project lbm_project)
