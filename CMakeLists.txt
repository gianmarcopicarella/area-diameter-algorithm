cmake_minimum_required (VERSION 3.18)

project(master-thesis)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

set(GOOGLETEST_PATH ${CMAKE_CURRENT_SOURCE_DIR}/external/googletest)

# dependencies
add_subdirectory(external/raylib)
add_subdirectory(external/googletest)
add_subdirectory(external/benchmark)

# main project
add_subdirectory(src)