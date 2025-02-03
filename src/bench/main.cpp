#include <benchmark/benchmark.h>

#include <vector>
#include <filesystem>
#include <iostream>
#include <chrono>
#include <type_traits>

#include "../common/Eppstein.h"
#include "../common/Antipodal.h"
#include "../common/Parser.h"

namespace fs = std::filesystem;

enum class Distribution
{
    UNIFORM,
    GAUSSIAN
};

enum class Algorithm
{
    EPPSTEIN,
    ANTIPODAL
};


static void Setup(const benchmark::State& state)
{
    static bool initCSV = false;
    if(!initCSV)
    {
        initCSV = true;
        constexpr auto shouldAppendToFile = false;
        MT::SZ::WriteLineToCsv("benchmark_data_results.csv", [](){
            return "name, count, area, diameter_indices, hull_indices";
        }, shouldAppendToFile);
    }
}


template<Distribution D>
fs::path GetPathToFolder(const benchmark::State& aState)
{
    const fs::path basePath {"../../data/samples/experiments"};
    if constexpr (D == Distribution::UNIFORM)
    {
        return basePath / fs::path{ "uniform/" + std::to_string(aState.range(0)) };
    }
    else
    {
        return basePath / fs::path{ "gaussian/" + std::to_string(aState.range(0)) };
    }
}

template<Distribution D, Algorithm A>
void BM_Template(benchmark::State& aState)
{
    constexpr auto maxAllowedArea = 4.l; // 2 mm^2
    constexpr auto maxAllowedDiameter = 2.l; // 4 mm
    constexpr auto maxAllowedPoints = (size_t) - 1; // No limit
    constexpr auto reconstructHull = true;

    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::time_point;

    size_t fileIndex = 0;
    for (auto _ : aState)
    {
        std::vector<MT::CM::Point2> points;
        const fs::path filename { "points_" + std::to_string(fileIndex) + ".csv" };
        MT::SZ::ReadPointsFromFile(GetPathToFolder<D>(aState) / filename, points);
        std::optional<MT::ConvexArea> resultOpt;

        auto start = std::chrono::high_resolution_clock::now();
        if constexpr (A == Algorithm::EPPSTEIN)
        {
            resultOpt = MT::EppsteinAlgorithm(points, maxAllowedPoints, maxAllowedArea, reconstructHull);
        }
        else
        {
            resultOpt = MT::AntipodalAlgorithm(points, maxAllowedPoints, maxAllowedArea, maxAllowedDiameter, reconstructHull);
        }

        auto end = std::chrono::high_resolution_clock::now();
        const auto elapsed =
                std::chrono::duration_cast<std::chrono::duration<double>>(end - start);

        MT::SZ::WriteLineToCsv("benchmark_data_results.csv", [&](){
            const auto& name = aState.name() + "/" + std::to_string(aState.range(0)) + "/iterations:" + std::to_string(fileIndex);
            MT::ConvexArea result = resultOpt.has_value() ? *resultOpt : MT::ConvexArea{};
            const auto& count = std::to_string(result.myPointsCount);
            const auto& area = std::to_string(result.myHullArea);
            const auto& diameter =
                    std::to_string(result.myDiameterOpt->myFirstIndex) + "\t" +
                    std::to_string(result.myDiameterOpt->mySecondIndex);

            std::string hullIndices = std::to_string(result.myHullIndices[0]);
            for(size_t i = 1; i < result.myHullIndices.size(); ++i)
            {
                hullIndices += "\t" + std::to_string(result.myHullIndices[i]);
            }
            return name + ", " + count + ", " + area + ", " + diameter + ", " + hullIndices;
        });

        ++fileIndex;
    }
}

BENCHMARK(BM_Template<Distribution::UNIFORM, Algorithm::EPPSTEIN>)
->Setup(Setup)->Name("Eppstein/Uniform")->Unit(benchmark::kMillisecond)->DenseRange(60, 70, 10)->Iterations(10);
/*
BENCHMARK(BM_Template<Distribution::GAUSSIAN, Algorithm::EPPSTEIN>)
->Setup(Setup)->Name("Eppstein/Gaussian")->Unit(benchmark::kMillisecond)->DenseRange(1, 10, 1)->Iterations(10);

BENCHMARK(BM_Template<Distribution::UNIFORM, Algorithm::ANTIPODAL>)
->Setup(Setup)->Name("Antipodal/Uniform")->Unit(benchmark::kMillisecond)->DenseRange(60, 150, 10)->Iterations(10);

BENCHMARK(BM_Template<Distribution::GAUSSIAN, Algorithm::ANTIPODAL>)
->Setup(Setup)->Name("Antipodal/Gaussian")->Unit(benchmark::kMillisecond)->DenseRange(1, 10, 1)->Iterations(10);
*/

BENCHMARK_MAIN();