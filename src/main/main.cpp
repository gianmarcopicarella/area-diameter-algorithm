#include <benchmark/benchmark.h>
#include <args.hxx>

#include <vector>
#include <nlohmann/json.hpp>

#include "../common/Eppstein.h"
#include "../common/Antipodal.h"
#include "../common/AntipodalOptimized.h"
#include "../common/Parser.h"
#include "../common/Constants.h"
#include "../common/TestingUtils.h"

#include <iostream>

int main(int argc, char** argv)
{
    constexpr auto defaultAlgorithm = static_cast<size_t>(Algorithm::EPPSTEIN);
    constexpr auto defaultMaxPointsCount = (size_t) - 1;
    constexpr auto defaultMaxArea = std::numeric_limits<long double>::infinity();
    constexpr auto defaultReconstructHull = false;

    args::ArgumentParser parser("CLI for Eppstein and Antipodal algorithms");
    args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});

    args::ValueFlag<std::string> filepath(parser,
                                          "Filepath",
                                          "Path to JSON file containing the point set",
                                          {'f', "Filepath"}, args::Options::Required);

    args::ValueFlag<size_t> algorithm(parser,
                                         "Algorithm",
                                         "Algorithm to be used: 0 = Eppstein, 1 = Antipodal, 2 = Antipodal Optimized",
                                         {'g', "Algorithm"}, defaultAlgorithm);

    args::ValueFlag<size_t> maxAllowedPoints(parser,
                                             "MaxCardinality",
                                             "Convex area's maximum allowed cardinality",
                                             {'m', "MaxCardinality"}, defaultMaxPointsCount);

    args::ValueFlag<long double> maxAllowedArea(parser,
                                                "MaxArea",
                                                "Convex area's maximum allowed area",
                                                {'a', "MaxArea"}, defaultMaxArea);

    args::Flag shouldReconstructHull(parser,
                                     "ReconstructHull",
                                     "If enabled, the algorithm returns a list of indices defining the optimal convex area",
                                     {'r', "ReconstructHull"}, defaultReconstructHull);

    args::ValueFlag<long double> maxAllowedDiameter(parser,
                                                    "MaxDiameter",
                                                    "Convex area's maximum allowed diameter", {'d', "MaxDiameter"});
    try
    {
        parser.ParseCLI(argc, argv);
    }
    catch (const args::Completion& e)
    {
        std::cout << e.what();
        return 0;
    }
    catch (const args::Help&)
    {
        std::cout << parser;
        return 0;
    }
    catch (const args::ParseError& e)
    {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }

    if(algorithm.Get() < 0 || algorithm.Get() >= static_cast<size_t>(Algorithm::COUNT))
    {
        std::cerr << "Error: invalid value for --Algorithm.\n";
        return 1;
    }
    else if(!maxAllowedDiameter && algorithm.Get() >= static_cast<size_t>(Algorithm::ANTIPODAL))
    {
        std::cerr << "Error: missing value for --MaxDiameter.\n";
        return 1;
    }
    else if(maxAllowedDiameter && algorithm.Get() < static_cast<size_t>(Algorithm::ANTIPODAL))
    {
        std::cerr << "Error: --MaxDiameter can only be specified when --Algorithm is set to 1 (Antipodal) or 2 (Antipodal Optimized).\n";
        return 1;
    }

    std::vector<MT::CM::Point2> points;
    MT::SZ::ReadPointsFromFile(filepath.Get(), points);
    std::optional<MT::ConvexArea> resultOpt;
    std::optional benchmarkInfo =  MT::BenchmarkInfo {};

    benchmark::ClobberMemory();
    StartHeapProfiling();

    const auto startTime = std::chrono::high_resolution_clock::now();
    benchmark::DoNotOptimize(points);

    switch (static_cast<Algorithm>(algorithm.Get()))
    {
        case Algorithm::EPPSTEIN:
        {
            resultOpt = MT::EppsteinAlgorithmWithBenchmarkInfo(points, benchmarkInfo, maxAllowedPoints.Get(), maxAllowedArea.Get(), shouldReconstructHull.Get());
            break;
        }
        case Algorithm::ANTIPODAL:
        {
            resultOpt = MT::AntipodalAlgorithmWithBenchmarkInfo(points, benchmarkInfo, maxAllowedPoints.Get(), maxAllowedArea.Get(), maxAllowedDiameter.Get(), shouldReconstructHull.Get());
            break;
        }
        case Algorithm::ANTIPODAL_OPTIMIZED:
        {
            resultOpt = MT::AntipodalOptimizedAlgorithmWithBenchmarkInfo(points, benchmarkInfo, maxAllowedPoints.Get(), maxAllowedArea.Get(), maxAllowedDiameter.Get(), shouldReconstructHull.Get());
            break;
        }
        default:
            break;
    }

    benchmark::DoNotOptimize(resultOpt);
    const auto endTime = std::chrono::high_resolution_clock::now();

    StopHeapProfiling();
    benchmark::ClobberMemory();

    json result;

    if(resultOpt)
    {
        MT::Solution solution;
        solution.myMaxCount = maxAllowedPoints.Get();
        solution.myMaxArea = maxAllowedArea.Get();
        if(algorithm.Get() >= static_cast<size_t>(Algorithm::ANTIPODAL))
        {
            solution.myMaxDiameter = maxAllowedDiameter.Get();
        }
        solution.myConvexAreaOpt = resultOpt;
        result["solution"] = solution;
    }

    result["time"] = std::chrono::duration<double, std::milli> {endTime - startTime}.count();
    result["memory"] = totalAllocatedBytes;
    result["entries"] = benchmarkInfo->myAllocatedEntriesCount;
    result["min_entries"] = benchmarkInfo->myRequiredEntriesCount;

    std::cout << result << std::endl;
    return 0;
}

