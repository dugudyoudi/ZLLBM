# /examples/lbm/CMakeLists.txt

if(BUILD_EXAMPLES)
    set(LIST_DIRECTORIES immersed_boundary)

    list(TRANSFORM LIST_DIRECTORIES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/")

    set(subdirs_with_cmake "")

    foreach(dir IN LISTS LIST_DIRECTORIES)
        message(STATUS "Build examples in: ${dir}")
        # Check if the subdirectory contains a CMakeLists.txt file
        file(GLOB_RECURSE cmake_files "${dir}/*/CMakeLists.txt")
        FOREACH(file_path ${cmake_files})
            GET_FILENAME_COMPONENT(dir_path ${file_path} PATH)
            add_subdirectory(${dir_path})
        ENDFOREACH()
        
    endforeach()
endif(BUILD_EXAMPLES)