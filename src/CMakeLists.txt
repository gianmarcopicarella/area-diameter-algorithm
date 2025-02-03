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

create_exec(main)
create_exec(bench)
create_exec(test)

set_target_properties(raylib gtest gtest_main benchmark benchmark_main uninstall PROPERTIES FOLDER HiddenTargets)

# target_include_directories(test PRIVATE ${FAST_CSV_HEADER_PATH})
# target_include_directories(bench PRIVATE ${FAST_CSV_HEADER_PATH})
# target_include_directories(main PRIVATE ${FAST_CSV_HEADER_PATH})

target_link_to_libs(main PRIVATE "raylib")
target_link_to_libs(bench PRIVATE "raylib;benchmark::benchmark")
target_link_to_libs(test PRIVATE "raylib;gtest_main")

pybind11_add_module(thesis ${CMAKE_CURRENT_SOURCE_DIR}/py-binding/main.cpp;${COMMON_ALL})
target_link_to_libs(thesis PRIVATE "raylib")