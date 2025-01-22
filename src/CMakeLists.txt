find_package(CGAL REQUIRED)

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

if ( CGAL_AUTO_LINK_ENABLED )
    target_link_to_exec(main PRIVATE "raylib;tsl::hopscotch_map")
    target_link_to_exec(bench PRIVATE "raylib;benchmark::benchmark;tsl::hopscotch_map")
    target_link_to_exec(test PRIVATE "raylib;gtest_main;tsl::hopscotch_map")
else()
    target_link_to_exec(main PRIVATE "raylib;CGAL::CGAL;CGAL::Data;tsl::hopscotch_map")
    target_link_to_exec(bench PRIVATE "raylib;benchmark::benchmark;CGAL::CGAL;CGAL::Data;tsl::hopscotch_map")
    target_link_to_exec(test PRIVATE "raylib;gtest_main;CGAL::CGAL;CGAL::Data;tsl::hopscotch_map")
endif()