#include <gtest/gtest.h>

#include "../common/Eppstein.h"
#include "../common/Antipodal.h"
#include "../common/Parser.h"

//#define VERBOSE

#include <iostream>
#include <iomanip>

#include <filesystem>
namespace fs = std::filesystem;

using namespace MT;

size_t CountFilesInDirectory(const std::string& aPath)
{
    size_t count = 0;
    for (const auto& entry : fs::directory_iterator(aPath))
    {
        ++count;
    }
    return count;
}

bool AreThereCollinearPoints(const std::vector<CM::Point2>& somePoints)
{
    for(const auto& pi : somePoints)
    {
        for (const auto& pj: somePoints)
        {
            for (const auto& pl: somePoints)
            {
                if (pi.index != pj.index && pi.index != pl.index && pj.index != pl.index)
                {
                    if(CM::AreCollinear(pi, pj, pl))
                    {
                        return true;
                    }
                }
            }
        }
    }
    return false;
}



TEST(EppsteinWithSolutions, BasicAssertions)
{
    const std::vector<std::string> paths = {
            "../../data/samples/eppstein/30",
            /*"../../data/samples/eppstein/40",
            "../../data/samples/eppstein/40_norm",
            "../../data/samples/eppstein/50",
            "../../data/samples/eppstein/80",
            "../../data/samples/eppstein/80_1"
             */
    };

    for (const auto & path : paths)
    {
        const auto filesCount = CountFilesInDirectory(path);
#ifdef VERBOSE
        std::cout << "TESTING SAMPLES IN FOLDER: " << path << std::endl;
#endif
        for(size_t i = 0; i < filesCount / 2; ++i)
        {
            const auto sampleFilename = std::string("data_") + std::to_string(i) + ".txt";
            const auto resultsFilename = std::string("results_") + std::to_string(i) + ".txt";
            std::vector<CM::Point2> points;
            std::vector<long double> solutionsPY;
            SZ::ReadPointsFromFile(fs::path{path} / fs::path{sampleFilename}, points);
            constexpr auto maxAllowedArea = std::numeric_limits<long double>::infinity();
            constexpr auto shouldReconstructHull = true;
            const auto res = MT::EppsteinAlgorithm(points, points.size(), maxAllowedArea, shouldReconstructHull);
            auto maxSolutionsDistance = 0.l;
            if(SZ::ReadSolutionsFromFile(fs::path{path} / fs::path{resultsFilename}, solutionsPY))
            {
                EXPECT_EQ(solutionsPY.size(), (res.results.size() - 3));
                for(int m = 0; m < solutionsPY.size(); ++m)
                {
                    const auto diff = res.results[m] - solutionsPY[m];
                    maxSolutionsDistance = std::max(maxSolutionsDistance, std::fabs(diff));
                    EXPECT_TRUE(CM::IsCloseToZero(diff, 0.001));
                }
            }
#ifdef VERBOSE
            const auto areThereCollinearPoints = AreThereCollinearPoints(points);
            std::cout   << std::fixed << std::setprecision(8) << "File " << (1+i) << "/" << (filesCount / 2) << ". Collinear points="
                        << (areThereCollinearPoints ? "yes" : "no") << ". Max distance=" << std::to_string(maxSolutionsDistance) << std::endl;
#endif
        }
    }
}


TEST(AntipodalWithSolutions, BasicAssertions)
{
    const std::vector<std::string> paths = {
            "../../data/samples/antipodal/60"
    };

    const auto sampleFilename = std::string("data_") + std::to_string(0) + ".txt";
    const auto resultsAreaFilename = std::string("results_areas_") + std::to_string(0) + ".txt";
    const auto resultsCountsFilename = std::string("results_counts_") + std::to_string(0) + ".txt";

    std::vector<CM::Point2> points;
    std::vector<long double> solutionsAreas;
    std::vector<size_t> solutionsCounts;

    SZ::ReadPointsFromFile(fs::path{paths[0]} / fs::path{sampleFilename}, points);
    SZ::ReadSolutionsFromFile(fs::path{paths[0]} / fs::path{resultsAreaFilename}, solutionsAreas);
    SZ::ReadSolutionsFromFile(fs::path{paths[0]} / fs::path{resultsCountsFilename}, solutionsCounts);

    constexpr auto maxAllowedArea = std::numeric_limits<long double>::infinity();
    constexpr auto maxAllowedDiameter = std::numeric_limits<long double>::infinity();
    const auto maxPointsCount = points.size();
    constexpr auto shouldReconstructHull = true;

    const auto res = MT::AntipodalAlgorithm(points, maxPointsCount, maxAllowedArea, maxAllowedDiameter, shouldReconstructHull);
    auto maxSolutionsDistance = 0.l;

    EXPECT_EQ(solutionsAreas.size(), res.results.size());
    EXPECT_EQ(solutionsCounts.size(), res.results.size());

    for(int m = 0; m < res.results.size(); ++m)
    {
        const auto diff = res.results[m].first - solutionsAreas[m];
        maxSolutionsDistance = std::max(maxSolutionsDistance, std::fabs(diff));
        std::cout << m << ", " << res.results[m].first << ", " << solutionsAreas[m] << ", " << res.results[m].second << ", " << solutionsCounts[m] << std::endl;
        EXPECT_TRUE(CM::IsCloseToZero(diff, 0.001) ||
                            (res.results[m].first == std::numeric_limits<long double>::infinity() &&
                             solutionsAreas[m] == std::numeric_limits<long double>::infinity()));
        EXPECT_EQ(res.results[m].second, solutionsCounts[m]);
    }
}

// Test with increasing k=3 to n with max diam = 50;
// Test with increasing k=3 to n with max diam = 100;
// Test with increasing k=3 to n with max diam = 400;

// Cross check Eppstein and Antipodal from k=3 to n. First run Eppstein and store the area and diameter. Then use the area and diameter as constraints in Antipodal.
// Area and Diameter should match
