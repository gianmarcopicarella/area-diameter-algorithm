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

function(target_link_to_exec exec visibility libs)
    if(TARGET ${exec})
        target_link_libraries(${exec} ${visibility} ${libs})
        message(STATUS "Added ${libs} to Executable ${exec}")
    endif() 
endfunction()

create_exec(main)
create_exec(bench)
create_exec(test)

set_target_properties(raylib gtest gtest_main benchmark benchmark_main uninstall PROPERTIES FOLDER HiddenTargets)

target_link_to_exec(main PRIVATE "raylib")
target_link_to_exec(bench PRIVATE "raylib;benchmark::benchmark")
target_link_to_exec(test PRIVATE "raylib;gtest_main")