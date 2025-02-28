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


# Path to folders
set(PATH_TO_REL_BUILD ${CMAKE_CURRENT_BINARY_DIR}) # Path to Release Build used in run_real_data_with_timeout.py
set(PATH_TO_DATA ${CMAKE_CURRENT_SOURCE_DIR}/../data) # Path to data folder
set(PATH_TO_EXPERIMENTS ${PATH_TO_DATA}/samples/experiments) # Path to experiments folder (contains uniform, gaussian and real samples)
set(PATH_TO_TESTS ${PATH_TO_DATA}/samples/tests) # Path to tests folder (used for unit tests in tests/main.cpp)
set(PATH_TO_REPORTS ${PATH_TO_DATA}/reports) # Path to reports folder
set(PATH_TO_PLOTS ${PATH_TO_DATA}/plots) # Path to plots folder (used to produce plots showing the solutions found during benchmarks)

# Path to files
set(PATH_TO_RAW_DATASET ${PATH_TO_DATA}/raw/detections_subset.json) # Raw dataset containing real-world data
set(PATH_TO_BENCHMARK_RUNS ${PATH_TO_REPORTS}/benchmark_runs.json) # Contains the data related to the program execution during benchmarks (time, memory, etc..)
set(PATH_TO_BENCHMARK_RESULTS ${PATH_TO_REPORTS}/benchmark_results.json) # Contains the data related to the solutions found during benchmarks with our solutions
set(PATH_TO_AREA_SELECTOR_RESULTS ${PATH_TO_REPORTS}/area_selector_results.json) # Contains the data related to the solutions found by the Area-Selector method
set(PATH_TO_LATEX_REPORT ${PATH_TO_DATA}/latex/main.tex) # Contains a LaTeX document presenting all the data collected during the experiments

# [Settings]
# Density and Standard deviation
set(DENSITY_MIN 100) # Starting number of points
set(DENSITY_STEP 10) # Increase amount at every step
set(DENSITIES_COUNT 11) # Number of experiments

set(STDDEV_MIN 0.5) # Starting standard deviation
set(STDDEV_STEP 0.5) # Increase amount at every step
set(STDDEVS_COUNT 13) # Number of experiments

# Experiments generation
set(SYNTHETIC_BENCHMARK_ITERATIONS 100) # Number of times the same experiments should be repeated
set(DETECTION_CONFIDENCE_THRESHOLD 0.84) # Detection probability threshold used to filter the raw dataset's detections
set(DETECTION_METHOD midog21_1st_stage) # Detection algorithm considered during the raw dataset processing
set(REAL_BENCHMARKS_COUNT 10) # Number of experiments using real world data (1 repetition by design)

# Diameters used for experiments
set(SYNTHETIC_BENCHMARK_DIAMETERS 2 3 4 5 6) # Diameters used for benchmarking with synthetic data
set(REAL_BENCHMARK_DIAMETERS 4.243) # Diameters used for benchmarking with real-world data
# additional processing required to use constexpr std::array<T, N> in C++
list(LENGTH SYNTHETIC_BENCHMARK_DIAMETERS SYNTHETIC_BENCHMARK_DIAMETERS_COUNT)
list(LENGTH REAL_BENCHMARK_DIAMETERS REAL_BENCHMARK_DIAMETERS_COUNT)
list(JOIN SYNTHETIC_BENCHMARK_DIAMETERS "," SYNTHETIC_BENCHMARK_DIAMETERS)
list(JOIN REAL_BENCHMARK_DIAMETERS "," REAL_BENCHMARK_DIAMETERS)


# Optimization
set(ENABLE_OPTIMIZATIONS_WITH_SYNTHETIC_DATA false) # Enable Antipodal and Eppstein's optimizations while benchmarking with synthetic data
set(ENABLE_OPTIMIZATIONS_WITH_REAL_DATA true) # Enable Antipodal and Eppstein's optimizations while benchmarking with real-world data

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/common/Constants.h.in ${CMAKE_CURRENT_SOURCE_DIR}/common/Constants.h @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/../python/constants.py.in ${CMAKE_CURRENT_SOURCE_DIR}/../python/constants.py @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/../python/constants.py.in ${CMAKE_CURRENT_SOURCE_DIR}/../external/Area-Selector/src/constants.py @ONLY)