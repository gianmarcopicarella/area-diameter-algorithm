#include <gtest/gtest.h>

#include "../common/Eppstein.h"
#include "../common/Antipodal.h"
#include "../common/Parser.h"

//#define VERBOSE

#include <iostream>
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
            "../../data/samples/eppstein/20",
            /*"../../data/samples/eppstein/30",
            "../../data/samples/eppstein/40",
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

            const auto& metadata = SZ::ReadSolutionsFromFile(fs::path{path} / fs::path{resultsFilename}, solutionsPY);

            EXPECT_TRUE(metadata.has_value());
            const auto maxAllowedArea = std::get<0>(metadata.value());
            /* const auto maxAllowedDiameter = std::get<1>(metadata.value()); */
            const auto maxAllowedPointsCount = std::get<2>(metadata.value());
            constexpr auto shouldReconstructHull = true;

            const auto res = MT::EppsteinAlgorithm(points, maxAllowedPointsCount, maxAllowedArea, shouldReconstructHull);

            auto maxSolutionsDistance = 0.l;
            EXPECT_EQ(solutionsPY.size(), (res.results.size() - 3));
            for(int m = 0; m < solutionsPY.size(); ++m)
            {
                const auto diff = res.results[m] - solutionsPY[m];
                maxSolutionsDistance = std::max(maxSolutionsDistance, std::fabs(diff));
                EXPECT_TRUE(CM::IsCloseToZero(diff, 0.001));
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
            SZ::ReadPointsFromFile(fs::path{path} / fs::path{sampleFilename}, points);

            std::vector<long double> correctAreas;
            std::vector<size_t> correctCounts;
            const auto& metadata = SZ::ReadSolutionsFromFile(fs::path{path} / fs::path{resultsFilename}, correctAreas, correctCounts);
            EXPECT_TRUE(metadata.has_value());

            const auto maxAllowedArea = std::get<0>(metadata.value());
            const auto maxAllowedDiameter = std::get<1>(metadata.value());
            const auto maxAllowedPointsCount = std::get<2>(metadata.value());
            constexpr auto shouldReconstructHull = true;

            const auto res = MT::AntipodalAlgorithm(points, maxAllowedPointsCount, maxAllowedArea, maxAllowedDiameter, shouldReconstructHull);

            auto maxSolutionsDistance = 0.l;
            EXPECT_EQ(correctAreas.size(), res.results.size());
            EXPECT_EQ(correctCounts.size(), res.results.size());

            for(int m = 0; m < res.results.size(); ++m)
            {
                const auto diff = res.results[m].first - correctAreas[m];
                maxSolutionsDistance = std::max(maxSolutionsDistance, std::fabs(diff));
                // std::cout << m << ", " << res.results[m].first << ", " << correctAreas[m] << ", " << res.results[m].second << ", " << correctCounts[m] << std::endl;
                EXPECT_TRUE(CM::IsCloseToZero(diff, 0.001) ||
                            (res.results[m].first == std::numeric_limits<long double>::infinity() &&
                                    correctAreas[m] == std::numeric_limits<long double>::infinity()));
                EXPECT_EQ(res.results[m].second, correctCounts[m]);
            }
#ifdef VERBOSE
            const auto areThereCollinearPoints = AreThereCollinearPoints(points);
            std::cout   << std::fixed << std::setprecision(8) << "File " << (1+i) << "/" << (filesCount / 2) << ". Collinear points="
                        << (areThereCollinearPoints ? "yes" : "no") << ". Max distance=" << std::to_string(maxSolutionsDistance) << std::endl;
#endif
        }
    }
}

// Test with increasing k=3 to n with max diam = 50;
// Test with increasing k=3 to n with max diam = 100;
// Test with increasing k=3 to n with max diam = 400;

// Cross check Eppstein and Antipodal from k=3 to n. First run Eppstein and store the area and diameter. Then use the area and diameter as constraints in Antipodal.
// Area and Diameter should match
