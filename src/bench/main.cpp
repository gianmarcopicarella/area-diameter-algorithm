#include <benchmark/benchmark.h>

#include <vector>
#include <filesystem>
#include <chrono>
#include <nlohmann/json.hpp>

#include "../common/Eppstein.h"
#include "../common/Antipodal.h"
#include "../common/Parser.h"
#include "../common/Constants.h"


using json = nlohmann::json;
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


template<Distribution D>
fs::path GetPathToFolder(const benchmark::State& aState)
{
    if constexpr (D == Distribution::UNIFORM)
    {
        return fs::path { MT::Constants::EXPERIMENT_SAMPLES_PATH } / fs::path{ "uniform/" + std::to_string(aState.range(0)) };
    }
    else
    {
        return fs::path { MT::Constants::EXPERIMENT_SAMPLES_PATH } / fs::path{ "gaussian/" + std::to_string(aState.range(0)) };
    }
}


template<Distribution D, Algorithm A>
void BM_Template(benchmark::State& aState, std::vector<MT::Solution>& someOutSolutions)
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
        const fs::path filename { "points_" + std::to_string(fileIndex) + ".json" };
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

        MT::Solution currentSolution;
        currentSolution.myId = aState.name() + "/" + std::to_string(aState.range(0)) + "/iterations:" + std::to_string(fileIndex);
        currentSolution.myMaxCount = maxAllowedPoints;
        currentSolution.myMaxArea = maxAllowedArea;
        currentSolution.myMaxDiameter = maxAllowedDiameter;
        currentSolution.myConvexAreaOpt = resultOpt;
        someOutSolutions.emplace_back(currentSolution);

        ++fileIndex;
    }
}


std::vector<MT::Solution> benchmarkSolutions;

BENCHMARK_TEMPLATE2_CAPTURE(BM_Template, Distribution::UNIFORM, Algorithm::EPPSTEIN, BM_Eppstein_Uniform, benchmarkSolutions)
->Name("Eppstein/Uniform")->Unit(benchmark::kMillisecond)->DenseRange(60, 150, 10)->Iterations(10);

BENCHMARK_TEMPLATE2_CAPTURE(BM_Template, Distribution::GAUSSIAN, Algorithm::EPPSTEIN, BM_Eppstein_Gaussian, benchmarkSolutions)
->Name("Eppstein/Gaussian")->Unit(benchmark::kMillisecond)->DenseRange(1, 10, 1)->Iterations(10);

BENCHMARK_TEMPLATE2_CAPTURE(BM_Template, Distribution::UNIFORM, Algorithm::ANTIPODAL, BM_Antipodal_Uniform, benchmarkSolutions)
->Name("Antipodal/Uniform")->Unit(benchmark::kMillisecond)->DenseRange(60, 150, 10)->Iterations(10);

BENCHMARK_TEMPLATE2_CAPTURE(BM_Template, Distribution::GAUSSIAN, Algorithm::ANTIPODAL, BM_Antipodal_Gaussian, benchmarkSolutions)
->Name("Antipodal/Gaussian")->Unit(benchmark::kMillisecond)->DenseRange(1, 10, 1)->Iterations(10);


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