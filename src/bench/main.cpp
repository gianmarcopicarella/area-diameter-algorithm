#include <benchmark/benchmark.h>

#include <vector>
#include <filesystem>
#include <nlohmann/json.hpp>

#include "../common/Eppstein.h"
#include "../common/Antipodal.h"
#include "../common/Parser.h"
#include "../common/Constants.h"
#include "../common/TestingUtils.h"
#include "../common/Constants.h"

#include <iostream>

void AddExtraCounter(const std::string& aName,
                     const std::vector<benchmark::IterationCount>& someValues,
                     benchmark::State& anOutState,
                     benchmark::Counter::Flags aFlags = benchmark::Counter::Flags::kDefaults,
                     benchmark::Counter::OneK aOneK = benchmark::Counter::OneK::kIs1000)
{
    const auto mean = Mean(someValues);
    anOutState.counters[aName + "_avg"] = benchmark::Counter(mean, aFlags, aOneK);
    anOutState.counters[aName + "_std"] = benchmark::Counter(StandardDeviation(mean, someValues), aFlags, aOneK);
}

namespace fs = std::filesystem;
template<Data D>
fs::path PathToData(const benchmark::State& aState)
{
    constexpr std::array<std::string_view, 3> folders = { "uniform", "gaussian", "real" };
    const std::string folder { folders[static_cast<size_t>(D)] };
    return fs::path { MT::Constants::PATH_TO_EXPERIMENTS } / fs::path{ folder } / fs::path { std::to_string(aState.range(0)) };
}

template<Algorithm A>
std::string ResultId(const benchmark::State& aState, const size_t aFileIndex)
{
    if constexpr(A == Algorithm::ANTIPODAL)
    {
        return aState.name() + "/" + std::to_string(aState.range(0)) + "/" + std::to_string(aState.range(1)) + "/iterations:" + std::to_string(aFileIndex);
    }
    else
    {
        return aState.name() + "/" + std::to_string(aState.range(0)) +  + "/iterations:" + std::to_string(aFileIndex);
    }
}

template<Data D, Algorithm A>
void BM_Template(benchmark::State& aState)
{
    constexpr auto maxAllowedArea = 4.l; // 2 mm^2
    constexpr auto maxAllowedPoints = (size_t) - 1; // No limit
    constexpr auto reconstructHull = true;

    size_t fileIndex = 0;

    std::vector<MT::Solution> solutions;
    solutions.reserve(aState.iterations());

    std::vector<benchmark::IterationCount> allocatedBytesCount;
    std::vector<benchmark::IterationCount> requiredEntriesCount;
    std::vector<benchmark::IterationCount> allocatedEntriesCount;

    allocatedBytesCount.reserve(aState.iterations());
    requiredEntriesCount.reserve(aState.iterations());
    allocatedEntriesCount.reserve(aState.iterations());

    for (auto _ : aState)
    {
        std::vector<MT::CM::Point2> points;
        const fs::path filename { "points_" + std::to_string(fileIndex) + ".json" };
        MT::SZ::ReadPointsFromFile(PathToData<D>(aState) / filename, points);
        std::optional<MT::ConvexArea> resultOpt;
        std::optional<MT::BenchmarkInfo> benchmarkInfo = MT::BenchmarkInfo {};

        benchmark::ClobberMemory();
        StartHeapProfiling();

        const auto startTime = std::chrono::high_resolution_clock::now();
        benchmark::DoNotOptimize(points);
        if constexpr (A == Algorithm::EPPSTEIN)
        {
            resultOpt = MT::EppsteinAlgorithmWithBenchmarkInfo(points, benchmarkInfo, maxAllowedPoints, maxAllowedArea, reconstructHull);
        }
        else
        {
            const auto maxAllowedDiameter = aState.range(1);
            resultOpt = MT::AntipodalAlgorithmWithBenchmarkInfo(points, benchmarkInfo, maxAllowedPoints, maxAllowedArea, maxAllowedDiameter, reconstructHull);
        }
        benchmark::DoNotOptimize(resultOpt);
        const auto endTime = std::chrono::high_resolution_clock::now();

        StopHeapProfiling();
        benchmark::ClobberMemory();

        allocatedBytesCount.emplace_back(totalAllocatedBytes);
        requiredEntriesCount.emplace_back(benchmarkInfo->myRequiredEntriesCount);
        allocatedEntriesCount.emplace_back(benchmarkInfo->myAllocatedEntriesCount);
        // std::cout << "[" << fileIndex << "] AB " << allocatedBytesCount[fileIndex] << ", RE " << requiredEntriesCount[fileIndex] << ", AE " << allocatedEntriesCount[fileIndex] << std::endl;

        aState.SetIterationTime(std::chrono::duration<double, std::milli> {endTime - startTime}.count());

        MT::Solution currentSolution;
        currentSolution.myId = ResultId<A>(aState, fileIndex);
        currentSolution.myMaxCount = maxAllowedPoints;
        currentSolution.myMaxArea = maxAllowedArea;

        if constexpr (A == Algorithm::ANTIPODAL)
        {
            currentSolution.myMaxDiameter = aState.range(1);
        }

        currentSolution.myConvexAreaOpt = resultOpt;
        solutions.emplace_back(currentSolution);
        ++fileIndex;
    }

    AddExtraCounter("mem", allocatedBytesCount, aState, benchmark::Counter::kDefaults, benchmark::Counter::kIs1024);
    AddExtraCounter("entries", allocatedEntriesCount, aState);
    AddExtraCounter("min_entries", requiredEntriesCount, aState);

    const auto& filepath = fs::path { MT::Constants::PATH_TO_BENCHMARK_CUSTOM_REPORT };
    std::ifstream oldBenchmarkSolutionsFile {filepath};
    json data;

    if(oldBenchmarkSolutionsFile.good())
    {
        data = json::parse(oldBenchmarkSolutionsFile);
        oldBenchmarkSolutionsFile.close();
    }
    else
    {
        data = json{{"results", std::vector<MT::Solution>{}}, {"count", 0}};
    }

    data["results"].insert(data["results"].end(), solutions);
    data["count"] = data["count"].get<size_t>() + solutions.size();

    std::ofstream newBenchmarkSolutionsFile(filepath, std::ios::trunc);
    assert(newBenchmarkSolutionsFile.is_open());
    if(!newBenchmarkSolutionsFile.is_open()) throw std::runtime_error("Cannot write benchmark results to JSON file!");
    newBenchmarkSolutionsFile << std::setw(4) << data;
    newBenchmarkSolutionsFile.close();
}

// 1) Uniform distribution, Increasing density [100, 200, step=10]

BENCHMARK_TEMPLATE2(BM_Template, Data::SYNTHETIC_UNIFORM, Algorithm::ANTIPODAL)
->Name("Antipodal/Uniform")->Unit(benchmark::kMillisecond)
->ArgsProduct({ benchmark::CreateDenseRange(0, 10, 1), benchmark::CreateDenseRange(2, 6, 1) })->Iterations(100);

BENCHMARK_TEMPLATE2(BM_Template, Data::SYNTHETIC_UNIFORM, Algorithm::EPPSTEIN)
->Name("Eppstein/Uniform")->Unit(benchmark::kMillisecond)->DenseRange(0, 10, 1)->Iterations(100);

// 2) Gaussian distribution, Increasing standard deviation [0.5, 5.5, step=0.5]

BENCHMARK_TEMPLATE2(BM_Template, Data::SYNTHETIC_GAUSSIAN, Algorithm::ANTIPODAL)
->Name("Antipodal/Gaussian")->Unit(benchmark::kMillisecond)
->ArgsProduct({ benchmark::CreateDenseRange(0, 10, 1), benchmark::CreateDenseRange(2, 6, 1) })->Iterations(100);

BENCHMARK_TEMPLATE2(BM_Template, Data::SYNTHETIC_GAUSSIAN, Algorithm::EPPSTEIN)
->Name("Eppstein/Gaussian")->Unit(benchmark::kMillisecond)->DenseRange(0, 10, 1)->Iterations(100);

// 3) Real world data [10 different samples] NO TIMEOUT
/*
BENCHMARK_TEMPLATE2(BM_Template, Data::REAL, Algorithm::ANTIPODAL)
->Name("Antipodal/Real")->Unit(benchmark::kMillisecond)
->ArgsProduct({ benchmark::CreateDenseRange(1, 9, 1), benchmark::CreateDenseRange(2, 2, 1) })->Iterations(1);
*/

int main(int argc, char** argv)
{
    char arg0_default[] = "benchmark";
    char* args_default = arg0_default;
    if (!argv)
    {
      argc = 1;
      argv = &args_default;
    }
    ::benchmark::Initialize(&argc, argv);
    if (::benchmark::ReportUnrecognizedArguments(argc, argv))
    {
        return 1;
    }
    fs::remove(fs::path { MT::Constants::PATH_TO_BENCHMARK_CUSTOM_REPORT });
    ::benchmark::RunSpecifiedBenchmarks();
    ::benchmark::Shutdown();

    return 0;
  }