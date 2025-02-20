file(   GLOB_RECURSE
        COMMON_ALL
        ${CMAKE_CURRENT_SOURCE_DIR}/common/*.h
        ${CMAKE_CURRENT_SOURCE_DIR}/common/*.cpp)

function(create_exec exec)
    file(   GLOB_RECURSE 
            ${exec}_ALL
            ${CMAKE_CURRENT_SOURCE_DIR}/${exec}/*.h
            ${CMAKE_CURRENT_SOURCE_DIR}/${exec}/*.cpp)
    list(LENGTH ${exec}_ALL LENGTH)
    if(LENGTH GREATER_EQUAL 1)
        add_executable(${exec} "${${exec}_ALL};${COMMON_ALL}")
    else()
        message(WARNING "Cannot create executable ${exec}. The list of files is empty!")
    endif() 
endfunction()

function(target_link_to_libs target visibility libs)
    if(TARGET ${target})
        target_link_libraries(${target} ${visibility} ${libs})
        message(STATUS "Added ${libs} to Executable ${target}")
    endif() 
endfunction()

function(add_target_coverage target)
    if(TARGET ${target} AND CMAKE_BUILD_TYPE STREQUAL "Debug")
        target_compile_options(${target} PRIVATE -coverage)
        target_link_options(${target} PRIVATE -coverage)
        message(STATUS "Enabled coverage for target ${target}")
    endif()
endfunction()

add_compile_options(-O3)

create_exec(main)
create_exec(bench)
create_exec(test)

set_target_properties(gtest gtest_main benchmark benchmark_main PROPERTIES FOLDER HiddenTargets)

target_link_to_libs(main PRIVATE "taywee::args;benchmark::benchmark;nlohmann_json::nlohmann_json")
target_link_to_libs(bench PRIVATE "benchmark::benchmark;nlohmann_json::nlohmann_json")
target_link_to_libs(test PRIVATE "gtest_main;nlohmann_json::nlohmann_json")

pybind11_add_module(thesis ${CMAKE_CURRENT_SOURCE_DIR}/py-binding/main.cpp;${COMMON_ALL})
target_link_to_libs(thesis PRIVATE "nlohmann_json::nlohmann_json")
set_target_properties(thesis PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../python/wrapper-build)

set(BENCHMARK_OUT_DATA_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../data/reports)
set(UNIT_TESTS_SAMPLES_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../data/samples/tests)
set(EXPERIMENT_SAMPLES_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../data/samples/experiments)
set(BENCHMARK_OUT_CUSTOM_FILENAME benchmark_data_results.json)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/common/Constants.h.in ${CMAKE_CURRENT_SOURCE_DIR}/common/Constants.h)