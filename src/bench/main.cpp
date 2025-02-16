#include <benchmark/benchmark.h>

#include <vector>
#include <filesystem>
#include <nlohmann/json.hpp>

#include "../common/Eppstein.h"
#include "../common/Antipodal.h"
#include "../common/Parser.h"
#include "../common/Constants.h"

#include <iostream>

enum class Data
{
    SYNTHETIC_UNIFORM = 0,
    SYNTHETIC_GAUSSIAN,
    REAL
};

enum class Algorithm
{
    EPPSTEIN = 0,
    ANTIPODAL
};

static benchmark::IterationCount totalAllocatedBytes { 0 };
static bool shouldTrackHeapMemory { false };

static void StartHeapProfiling()
{
    totalAllocatedBytes = 0;
    shouldTrackHeapMemory = true;
}

static void StopHeapProfiling()
{
    shouldTrackHeapMemory = false;
}

void* operator new(size_t sz)
{
    if(shouldTrackHeapMemory)
    {
        totalAllocatedBytes += sz;
    }
    return std::malloc(sz);
}

template <typename T>
T Mean(const std::vector<T>& someValues)
{
    const auto sum = std::accumulate(someValues.begin(), someValues.end(), T { 0 });
    return sum / someValues.size();
}

template <typename T>
T StandardDeviation(const T& aMean, const std::vector<T>& someValues)
{
    T sumOfSquaredResiduals { 0 };
    for(const auto& aValue : someValues)
    {
        sumOfSquaredResiduals += (aMean - aValue) * (aMean - aValue);
    }
    return std::sqrtl(sumOfSquaredResiduals / someValues.size());
}

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
    constexpr std::array<std::string_view, 3> folders = { "uniform/", "gaussian/", "real/" };
    const std::string folder { folders[static_cast<size_t>(D)] };
    return fs::path { MT::Constants::EXPERIMENT_SAMPLES_PATH } / fs::path{ folder + std::to_string(aState.range(0)) };
}

template<Data D, Algorithm A>
void BM_Template(benchmark::State& aState, std::vector<MT::Solution>& someOutSolutions)
{
    constexpr auto maxAllowedArea = 4.l; // 2 mm^2
    constexpr auto maxAllowedPoints = (size_t) - 1; // No limit
    constexpr auto reconstructHull = true;

    size_t fileIndex = 0;

    std::vector<benchmark::IterationCount> allocatedBytesCount { aState.iterations(), 0 };
    std::vector<benchmark::IterationCount> requiredEntriesCount { aState.iterations(), 0 };
    std::vector<benchmark::IterationCount> allocatedEntriesCount { aState.iterations(), 0 };

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

        allocatedBytesCount[fileIndex] = totalAllocatedBytes;
        requiredEntriesCount[fileIndex] = benchmarkInfo->myRequiredEntriesCount;
        allocatedEntriesCount[fileIndex] = benchmarkInfo->myAllocatedEntriesCount;
        std::cout << "AB " << allocatedBytesCount[fileIndex] << ", RE " << requiredEntriesCount[fileIndex] << ", AE " << allocatedEntriesCount[fileIndex] << std::endl;

        aState.SetIterationTime(std::chrono::duration<double, std::milli> {endTime - startTime}.count());

        MT::Solution currentSolution;
        currentSolution.myId = aState.name() + "/" + std::to_string(aState.range(0)) + "/iterations:" + std::to_string(fileIndex);
        currentSolution.myMaxCount = maxAllowedPoints;
        currentSolution.myMaxArea = maxAllowedArea;

        if constexpr (A == Algorithm::ANTIPODAL)
        {
            currentSolution.myMaxDiameter = aState.range(1);
        }

        currentSolution.myConvexAreaOpt = resultOpt;
        someOutSolutions.emplace_back(currentSolution);
        ++fileIndex;
    }

    AddExtraCounter("mem", allocatedBytesCount, aState, benchmark::Counter::kDefaults, benchmark::Counter::kIs1024);
    AddExtraCounter("entries", allocatedEntriesCount, aState);
    AddExtraCounter("min_entries", requiredEntriesCount, aState);
}

std::vector<MT::Solution> benchmarkSolutions;

// 1) Uniform distribution, Increasing density [100, 200, step=10]

BENCHMARK_TEMPLATE2_CAPTURE(BM_Template, Data::SYNTHETIC_UNIFORM, Algorithm::ANTIPODAL, BM_Antipodal_Uniform, benchmarkSolutions)
->Name("Antipodal/Uniform")->Unit(benchmark::kMillisecond)->DenseRange(0, 10, 1)->DenseRange(2, 5, 1)->Iterations(100);
BENCHMARK_TEMPLATE2_CAPTURE(BM_Template, Data::SYNTHETIC_UNIFORM, Algorithm::EPPSTEIN, BM_Eppstein_Uniform, benchmarkSolutions)
->Name("Eppstein/Uniform")->Unit(benchmark::kMillisecond)->DenseRange(0, 10, 1)->Iterations(100);


// 2) Gaussian distribution, Increasing standard deviation [0.5, 3, step=0.25]
/*
BENCHMARK_TEMPLATE2_CAPTURE(BM_Template, Data::SYNTHETIC_GAUSSIAN, Algorithm::ANTIPODAL, BM_Antipodal_Gaussian, benchmarkSolutions)
->Name("Antipodal/Gaussian")->Unit(benchmark::kMillisecond)->DenseRange(0, 10, 1)->DenseRange(2, 5, 1)->Iterations(100);
BENCHMARK_TEMPLATE2_CAPTURE(BM_Template, Data::SYNTHETIC_GAUSSIAN, Algorithm::EPPSTEIN, BM_Eppstein_Gaussian, benchmarkSolutions)
->Name("Eppstein/Gaussian")->Unit(benchmark::kMillisecond)->DenseRange(0, 10, 1)->Iterations(100);
*/

// 3) Real world data [10 different samples]
/*
BENCHMARK_TEMPLATE2_CAPTURE(BM_Template, Data::REAL, Algorithm::ANTIPODAL, BM_Antipodal_Real, benchmarkSolutions)
->Name("Antipodal/Real")->Unit(benchmark::kMillisecond)->DenseRange(1, 9, 1)->DenseRange(2, 5, 1)->Iterations(1);
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
    ::benchmark::RunSpecifiedBenchmarks();
    ::benchmark::Shutdown();

    const auto& filepath = fs::path { MT::Constants::BENCHMARK_OUT_DATA_PATH } / fs::path{ "benchmark_data_results.json" };
    std::ofstream file(filepath, std::ios::trunc);
    assert(file.is_open());
    if(!file.is_open()) throw std::runtime_error("Cannot write benchmark results to JSON file!");
    file << std::setw(4) << json{{"results", benchmarkSolutions}, {"count", benchmarkSolutions.size()}};
    file.close();

    return 0;
  }