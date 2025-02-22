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
        constexpr std::string_view PATH_TO_BENCHMARK_CUSTOM_REPORT = "/Users/gianmarcopicarella/master-thesis/src/../data/reports/benchmark_data_results.json";
        constexpr std::string_view PATH_TO_BENCHMARK_BASE_REPORT = "/Users/gianmarcopicarella/master-thesis/src/../data/reports/benchmark_runs.json";
    }
}

#endif //MASTER_THESIS_CONSTANTS_H_IN_H
