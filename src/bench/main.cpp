#include <benchmark/benchmark.h>

#include <vector>
#include <filesystem>
#include <nlohmann/json.hpp>

#include "../common/Eppstein.h"
#include "../common/Antipodal.h"
#include "../common/Parser.h"
#include "../common/Constants.h"

#include <iostream>

#define MEMORY_PROFILER

#ifdef MEMORY_PROFILER
static benchmark::IterationCount totalAllocatedBytes { 0 };
static benchmark::IterationCount totalDeallocatedBytes { 0 };
static benchmark::IterationCount maxBytesUsed { 0 };
static bool shouldTrackHeapMemory { false };

static void StartHeapProfiling()
{
    totalAllocatedBytes = 0;
    totalDeallocatedBytes = 0;
    maxBytesUsed = 0;
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
        maxBytesUsed = std::max(maxBytesUsed, totalAllocatedBytes - totalDeallocatedBytes);
    }
    auto ptr = static_cast<size_t*>(std::malloc(sz+1));
    *ptr = sz;
    return ptr + 1;
}

void operator delete(void* ptr) noexcept
{
    auto b_ptr = static_cast<size_t*>(ptr) - 1;
    if(shouldTrackHeapMemory)
    {
        totalDeallocatedBytes += *b_ptr;
    }
    std::free(b_ptr);
}
#endif

using json = nlohmann::json;
namespace fs = std::filesystem;

enum class Distribution
{
    UNIFORM = 0,
    GAUSSIAN
};

enum class Algorithm
{
    EPPSTEIN = 0,
    ANTIPODAL
};

template <typename T>
T ComputeMean(const std::vector<T>& someValues)
{
    const auto sum = std::accumulate(someValues.begin(), someValues.end(), T { 0 });
    return sum / someValues.size();
}

template <typename T>
T ComputeSTDev(const T& aMean, const std::vector<T>& someValues)
{
    T sumOfSquaredResiduals { 0 };
    for(const auto& aValue : someValues)
    {
        sumOfSquaredResiduals += (aMean - aValue) * (aMean - aValue);
    }
    return std::sqrtl(sumOfSquaredResiduals / someValues.size());
}

void ComputeExtraMetrics(
        const std::vector<benchmark::IterationCount>& someMaxBytes,
        const std::vector<benchmark::IterationCount>& someMaxEntries,
        const std::vector<benchmark::IterationCount>& someMinEntries,
        benchmark::State& anOutState)
{
    if(!someMaxBytes.empty())
    {
        const auto maxMemMean = ComputeMean(someMaxBytes);
        const auto maxMemStd = ComputeSTDev(maxMemMean, someMaxBytes);
        anOutState.counters["avg_mem"] = benchmark::Counter(maxMemMean,
                                                            benchmark::Counter::kDefaults, benchmark::Counter::kIs1024);
        anOutState.counters["std_mem"] = benchmark::Counter(maxMemStd,
                                                            benchmark::Counter::kDefaults, benchmark::Counter::kIs1024);
    }
    else
    {
        anOutState.counters["avg_mem"] = -1;
        anOutState.counters["std_mem"] = -1;
    }

    if(!someMaxEntries.empty())
    {
        const auto entriesCountMean = ComputeMean(someMaxEntries);
        const auto entriesCountStd = ComputeSTDev(entriesCountMean, someMaxEntries);
        anOutState.counters["avg_entries"] = benchmark::Counter(entriesCountMean);
        anOutState.counters["std_entries"] = benchmark::Counter(entriesCountStd);
    }
    else
    {
        anOutState.counters["avg_entries"] = -1;
        anOutState.counters["std_entries"] = -1;
    }

    if(!someMinEntries.empty())
    {
        const auto entriesCountMean = ComputeMean(someMinEntries);
        const auto entriesCountStd = ComputeSTDev(entriesCountMean, someMinEntries);
        anOutState.counters["avg_min_entries"] = benchmark::Counter(entriesCountMean);
        anOutState.counters["std_min_entries"] = benchmark::Counter(entriesCountStd);
    }
    else
    {
        anOutState.counters["avg_min_entries"] = -1;
        anOutState.counters["std_min_entries"] = -1;
    }
}

template<Distribution D>
fs::path GetPathToExperiments(const benchmark::State& aState)
{
    constexpr std::array<std::string_view, 2> folders = { "uniform/", "gaussian/" };
    const std::string folder { folders[static_cast<size_t>(D)] };
    return fs::path { MT::Constants::EXPERIMENT_SAMPLES_PATH } / fs::path{ folder + std::to_string(aState.range(0)) };
}

template<Distribution D, Algorithm A>
void BM_Template(benchmark::State& aState, std::vector<MT::Solution>& someOutSolutions)
{
    constexpr auto maxAllowedArea = 4.l; // 2 mm^2
    constexpr auto maxAllowedDiameter = 2.l; // 4 mm
    constexpr auto maxAllowedPoints = (size_t) - 1; // No limit
    constexpr auto reconstructHull = true;

    size_t fileIndex = 0;
#ifdef MEMORY_PROFILER
    std::vector<benchmark::IterationCount> maxBytesUsedForAllRepetitions;
    std::vector<benchmark::IterationCount> maxEntriesUsedForAllRepetitions;
    std::vector<benchmark::IterationCount> minEntriesNeededForAllRepetitions;
    maxBytesUsedForAllRepetitions.reserve(aState.iterations());
    maxEntriesUsedForAllRepetitions.reserve(aState.iterations());
    minEntriesNeededForAllRepetitions.reserve(aState.iterations());
#endif
    for (auto _ : aState)
    {
        std::vector<MT::CM::Point2> points;
        const fs::path filename { "points_" + std::to_string(fileIndex) + ".json" };
        MT::SZ::ReadPointsFromFile(GetPathToExperiments<D>(aState) / filename, points);
        std::optional<MT::ConvexArea> resultOpt;
        std::optional<MT::BenchmarkInfo> benchmarkInfo = MT::BenchmarkInfo {};

#ifdef MEMORY_PROFILER
#if defined(__clang__)
        asm volatile("" ::: "memory");
#endif
        StartHeapProfiling();
#endif
        const auto startTime = std::chrono::high_resolution_clock::now();
        benchmark::DoNotOptimize(points);
        if constexpr (A == Algorithm::EPPSTEIN)
        {
            resultOpt = MT::EppsteinAlgorithmWithBenchmarkInfo(points, benchmarkInfo, maxAllowedPoints, maxAllowedArea, reconstructHull);
        }
        else
        {
            resultOpt = MT::AntipodalAlgorithmWithBenchmarkInfo(points, benchmarkInfo, maxAllowedPoints, maxAllowedArea, maxAllowedDiameter, reconstructHull);
        }
        benchmark::DoNotOptimize(resultOpt);
        const auto endTime = std::chrono::high_resolution_clock::now();

#ifdef MEMORY_PROFILER
        StopHeapProfiling();
#if defined(__clang__)
        asm volatile("" ::: "memory");
#endif
        if(benchmarkInfo->myMinEntriesCount > -1)
        {
            minEntriesNeededForAllRepetitions.emplace_back(benchmarkInfo->myMinEntriesCount);
        }
        if(benchmarkInfo->myCreatedEntriesCount > -1)
        {
            maxEntriesUsedForAllRepetitions.emplace_back(benchmarkInfo->myCreatedEntriesCount);
        }
        maxBytesUsedForAllRepetitions.emplace_back(maxBytesUsed);
        std::cout << "M " << maxBytesUsed << ", A " << totalAllocatedBytes << ", D " << totalDeallocatedBytes << ", S " << (totalAllocatedBytes - totalDeallocatedBytes) << std::endl;
#endif
        aState.SetIterationTime(std::chrono::duration<double, std::milli> {endTime - startTime}.count());

        MT::Solution currentSolution;
        currentSolution.myId = aState.name() + "/" + std::to_string(aState.range(0)) + "/iterations:" + std::to_string(fileIndex);
        currentSolution.myMaxCount = maxAllowedPoints;
        currentSolution.myMaxArea = maxAllowedArea;
        currentSolution.myMaxDiameter = maxAllowedDiameter;
        currentSolution.myConvexAreaOpt = resultOpt;
        someOutSolutions.emplace_back(currentSolution);
        ++fileIndex;
    }
#ifdef MEMORY_PROFILER
    ComputeExtraMetrics(maxBytesUsedForAllRepetitions, maxEntriesUsedForAllRepetitions, minEntriesNeededForAllRepetitions, aState);
#endif
}


std::vector<MT::Solution> benchmarkSolutions;

// 1) Uniform distribution, Increasing density [150, 250, step=10]
BENCHMARK_TEMPLATE2_CAPTURE(BM_Template, Distribution::UNIFORM, Algorithm::ANTIPODAL, BM_Antipodal_Uniform, benchmarkSolutions)
->Name("Antipodal/Uniform")->Unit(benchmark::kMillisecond)->DenseRange(0, 10, 1)->Iterations(10);
BENCHMARK_TEMPLATE2_CAPTURE(BM_Template, Distribution::UNIFORM, Algorithm::EPPSTEIN, BM_Eppstein_Uniform, benchmarkSolutions)
->Name("Eppstein/Uniform")->Unit(benchmark::kMillisecond)->DenseRange(0, 10, 1)->Iterations(10);

// 2) Gaussian distribution, Increasing standard deviation [0.5, 3, step=0.5]
BENCHMARK_TEMPLATE2_CAPTURE(BM_Template, Distribution::GAUSSIAN, Algorithm::ANTIPODAL, BM_Antipodal_Gaussian, benchmarkSolutions)
->Name("Antipodal/Gaussian")->Unit(benchmark::kMillisecond)->DenseRange(0, 10, 1)->Iterations(10);
BENCHMARK_TEMPLATE2_CAPTURE(BM_Template, Distribution::GAUSSIAN, Algorithm::EPPSTEIN, BM_Eppstein_Gaussian, benchmarkSolutions)
->Name("Eppstein/Gaussian")->Unit(benchmark::kMillisecond)->DenseRange(0, 10, 1)->Iterations(10);


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