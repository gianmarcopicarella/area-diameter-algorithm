//
// Created by Gianmarco Picarella on 03/02/25.
//

#ifndef MASTER_THESIS_CONSTANTS_H_IN_H
#define MASTER_THESIS_CONSTANTS_H_IN_H

#include <string_view>

namespace MT
{
    namespace Constants
    {
        constexpr std::string_view PATH_TO_TESTS = "/Users/gianmarcopicarella/master-thesis/src/../data/samples/tests";
        constexpr std::string_view PATH_TO_EXPERIMENTS = "/Users/gianmarcopicarella/master-thesis/src/../data/samples/experiments";
        constexpr std::string_view PATH_TO_BENCHMARK_RUNS = "/Users/gianmarcopicarella/master-thesis/src/../data/reports/benchmark_runs.json";
        constexpr std::string_view PATH_TO_BENCHMARK_RESULTS = "/Users/gianmarcopicarella/master-thesis/src/../data/reports/benchmark_results.json";

        constexpr std::array<long double, 5> SYNTHETIC_BENCHMARK_DIAMETERS = {2,3,4,5,6};

        constexpr size_t SYNTHETIC_BENCHMARK_UNIFORM_ITERATIONS = 100;
        constexpr size_t SYNTHETIC_BENCHMARK_GAUSSIAN_ITERATIONS = 100;

        constexpr size_t DENSITIES_COUNT = 11;
        constexpr size_t STDDEVS_COUNT = 13;

        constexpr bool ENABLE_OPTIMIZATIONS_WITH_SYNTHETIC_DATA = false;
    }
}

#endif //MASTER_THESIS_CONSTANTS_H_IN_H
