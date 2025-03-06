cmake_minimum_required (VERSION 3.18)

project(master-thesis)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

set(BENCHMARK_ENABLE_TESTING OFF CACHE INTERNAL "")
set(ARGS_BUILD_EXAMPLE OFF CACHE INTERNAL "")
set(ARGS_BUILD_UNITTESTS OFF CACHE INTERNAL "")
set(JSON_BuildTests OFF CACHE INTERNAL "")
set(PYTHON_EXECUTABLE /usr/bin/python3)

# dependencies
add_subdirectory(external/googletest)
add_subdirectory(external/benchmark)
add_subdirectory(external/pybind11)
add_subdirectory(external/json)
add_subdirectory(external/args)

# main project
add_subdirectory(src)