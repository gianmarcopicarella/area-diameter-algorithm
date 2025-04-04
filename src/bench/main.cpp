#include <benchmark/benchmark.h>

#include <vector>
#include <filesystem>
#include <nlohmann/json.hpp>

#include "../common/Eppstein.h"
#include "../common/Antipodal.h"
#include "../common/AntipodalOptimized.h"
#include "../common/Parser.h"
#include "../common/Constants.h"
#include "../common/TestingUtils.h"
#include "../common/Constants.h"

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
std::string ResultName(const benchmark::State& aState, const size_t aFileIndex)
{
    if constexpr(A >= Algorithm::ANTIPODAL)
    {
        return aState.name() + "/" + std::to_string(aState.range(0)) + "/" + std::to_string(aState.range(1)) + "/iterations:" + std::to_string(aFileIndex);
    }
    else
    {
        return aState.name() + "/" + std::to_string(aState.range(0)) +  + "/iterations:" + std::to_string(aFileIndex);
    }
}

template<Data D>
long double GetMaxDiameter(const benchmark::State& aState)
{
    return MT::Constants::SYNTHETIC_BENCHMARK_DIAMETERS[aState.range(1)];
}

template<Data D, Algorithm A>
void BM_Template(benchmark::State& aState, bool aShouldUseOptimizations)
{
    constexpr auto maxAllowedArea = 4.l; // 4 mm^2
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
        std::optional benchmarkInfo = MT::BenchmarkInfo {};

        benchmark::ClobberMemory();
        StartHeapProfiling();

        const auto startTime = std::chrono::high_resolution_clock::now();
        benchmark::DoNotOptimize(points);
        if constexpr (A == Algorithm::EPPSTEIN)
        {
            resultOpt = MT::EppsteinAlgorithmWithBenchmarkInfo(points, benchmarkInfo, maxAllowedPoints, maxAllowedArea, reconstructHull, aShouldUseOptimizations);
        }
        else
        {
            const auto maxAllowedDiameter = GetMaxDiameter<D>(aState);
            if constexpr(A == Algorithm::ANTIPODAL)
            {
                resultOpt = MT::AntipodalAlgorithmWithBenchmarkInfo(points, benchmarkInfo, maxAllowedPoints, maxAllowedArea, maxAllowedDiameter, reconstructHull, aShouldUseOptimizations);
            }
            else if constexpr(A == Algorithm::ANTIPODAL_OPTIMIZED)
            {
                resultOpt = MT::AntipodalOptimizedAlgorithmWithBenchmarkInfo(points, benchmarkInfo, maxAllowedPoints, maxAllowedArea, maxAllowedDiameter, reconstructHull, aShouldUseOptimizations);
            }
        }
        benchmark::DoNotOptimize(resultOpt);
        const auto endTime = std::chrono::high_resolution_clock::now();

        StopHeapProfiling();
        benchmark::ClobberMemory();

        allocatedBytesCount.emplace_back(totalAllocatedBytes);
        requiredEntriesCount.emplace_back(benchmarkInfo->myRequiredEntriesCount);
        allocatedEntriesCount.emplace_back(benchmarkInfo->myAllocatedEntriesCount);

        aState.SetIterationTime(std::chrono::duration<double, std::milli> {endTime - startTime}.count());

        MT::Solution currentSolution;
        currentSolution.myName = ResultName<A>(aState, fileIndex);
        currentSolution.myMaxCount = maxAllowedPoints;
        currentSolution.myMaxArea = maxAllowedArea;

        if constexpr (A >= Algorithm::ANTIPODAL)
        {
            currentSolution.myMaxDiameter = GetMaxDiameter<D>(aState);
        }

        currentSolution.myConvexAreaOpt = resultOpt;
        solutions.emplace_back(currentSolution);
        ++fileIndex;
    }

    AddExtraCounter("mem", allocatedBytesCount, aState, benchmark::Counter::kDefaults, benchmark::Counter::kIs1024);
    AddExtraCounter("entries", allocatedEntriesCount, aState);
    AddExtraCounter("min_entries", requiredEntriesCount, aState);

    const auto& filepath = fs::path { MT::Constants::PATH_TO_BENCHMARK_RESULTS };
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
BENCHMARK_TEMPLATE2_CAPTURE(BM_Template, Data::SYNTHETIC_UNIFORM, Algorithm::ANTIPODAL_OPTIMIZED, Antipodal_Uniform, MT::Constants::ENABLE_OPTIMIZATIONS_WITH_SYNTHETIC_DATA)
->Name("Antipodal/Uniform")->Unit(benchmark::kMillisecond)
->ArgsProduct({ benchmark::CreateDenseRange(0, MT::Constants::DENSITIES_COUNT - 1, 1),
                benchmark::CreateDenseRange(0, MT::Constants::SYNTHETIC_BENCHMARK_DIAMETERS.size() - 1, 1) })
->Iterations(MT::Constants::SYNTHETIC_BENCHMARK_ITERATIONS);

BENCHMARK_TEMPLATE2_CAPTURE(BM_Template, Data::SYNTHETIC_UNIFORM, Algorithm::EPPSTEIN, Eppstein_Uniform, MT::Constants::ENABLE_OPTIMIZATIONS_WITH_SYNTHETIC_DATA)
->Name("Eppstein/Uniform")->Unit(benchmark::kMillisecond)->DenseRange(0, MT::Constants::DENSITIES_COUNT - 1, 1)
->Iterations(MT::Constants::SYNTHETIC_BENCHMARK_ITERATIONS);


// 2) Gaussian distribution, Increasing standard deviation [0.5, 6.5, step=0.5]
BENCHMARK_TEMPLATE2_CAPTURE(BM_Template, Data::SYNTHETIC_GAUSSIAN, Algorithm::ANTIPODAL_OPTIMIZED, Antipodal_Gaussian, MT::Constants::ENABLE_OPTIMIZATIONS_WITH_SYNTHETIC_DATA)
->Name("Antipodal/Gaussian")->Unit(benchmark::kMillisecond)
->ArgsProduct({ benchmark::CreateDenseRange(0, MT::Constants::STDDEVS_COUNT - 1, 1),
                benchmark::CreateDenseRange(0, MT::Constants::SYNTHETIC_BENCHMARK_DIAMETERS.size() - 1, 1) })
->Iterations(MT::Constants::SYNTHETIC_BENCHMARK_ITERATIONS);

BENCHMARK_TEMPLATE2_CAPTURE(BM_Template, Data::SYNTHETIC_GAUSSIAN, Algorithm::EPPSTEIN, Eppstein_Gaussian, MT::Constants::ENABLE_OPTIMIZATIONS_WITH_SYNTHETIC_DATA)
->Name("Eppstein/Gaussian")->Unit(benchmark::kMillisecond)->DenseRange(0, MT::Constants::STDDEVS_COUNT - 1, 1)
->Iterations(MT::Constants::SYNTHETIC_BENCHMARK_ITERATIONS);


int main(int argc, char** argv)
{
    const std::vector<std::string> customArgs = {
            "--benchmark_out_format=json",
            "--benchmark_counters_tabular=true",
            "--benchmark_out=" + std::string{ MT::Constants::PATH_TO_BENCHMARK_RUNS }
    };
    if (argc == 1)
    {
        argc = customArgs.size();
        argv = static_cast<char **>(malloc(argc * sizeof(char *)));
        for(size_t i = 0; i < argc; ++i)
        {
            argv[i] = const_cast<char*>(customArgs[i].c_str());
        }
    }
    ::benchmark::Initialize(&argc, argv);
    if (::benchmark::ReportUnrecognizedArguments(argc, argv))
    {
        return 1;
    }
    fs::remove(fs::path { MT::Constants::PATH_TO_BENCHMARK_RESULTS });
    ::benchmark::RunSpecifiedBenchmarks();
    ::benchmark::Shutdown();

    return 0;
  }