#include <gtest/gtest.h>

#include "../common/Eppstein.h"
#include "../common/Antipodal.h"
#include "../common/Parser.h"
#include "../common/Utils.h"

#define VERBOSE

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
                if (pi.myIndex != pj.myIndex && pi.myIndex != pl.myIndex && pj.myIndex != pl.myIndex)
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

void GetHullPoints(const std::vector<size_t>& someIndices,
                   const std::vector<CM::Point2>& somePoints,
                   std::vector<CM::Point2>& someOutHullPoints)
{
    someOutHullPoints.resize(someIndices.size());
    std::transform(someIndices.begin(), someIndices.end(), someOutHullPoints.begin(),
                   [&](auto index){ return somePoints[index]; });
}

template<typename T>
void PrintResults(const T& aResult, const Diameter& aDiameter, const std::string& aName, bool anAddEndOfLineFlag = false)
{
    std::cout << aName;
    std::cout <<    " (Area= "   << aResult.myHullArea <<
                    ", Diam= "   << std::sqrtl(aDiameter.Value2()) <<
                    ", Pts= "    << aResult.myHullIndices.size() << "/" << aResult.myPointsCount << ")";
    if(anAddEndOfLineFlag)
    {
        std::cout << std::endl;
    }
}

TEST(EppsteinWithSolutions, BasicAssertions)
{
    const std::vector<std::string> paths = {
            "../../data/samples/eppstein/50",
            "../../data/samples/eppstein/60",
            "../../data/samples/eppstein/70",
            "../../data/samples/eppstein/80",
            "../../data/samples/eppstein/90",
            "../../data/samples/eppstein/100",
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
            std::vector<long double> correctAreas;
            SZ::ReadPointsFromFile(fs::path{path} / fs::path{sampleFilename}, points);

            const auto& metadata = SZ::ReadSolutionsFromFile(fs::path{path} / fs::path{resultsFilename}, correctAreas);

            EXPECT_TRUE(metadata.has_value());
            const auto maxAllowedArea = std::get<0>(metadata.value());
            /* const auto maxAllowedDiameter = std::get<1>(metadata.value()); */
            const auto maxAllowedPointsCount = std::get<2>(metadata.value());
            constexpr auto shouldReconstructHull = true;

            const auto res = MT::EppsteinAlgorithm(points, maxAllowedPointsCount, maxAllowedArea, shouldReconstructHull);

            auto maxSolutionsDistance = 0.l;
            EXPECT_EQ(correctAreas.size(), res.results.size());
            for(int m = 0; m < res.results.size(); ++m)
            {
                const auto diff = res.results[m] - correctAreas[m];
                EXPECT_TRUE((res.results[m] == std::numeric_limits<long double>::infinity() &&
                            correctAreas[m] == std::numeric_limits<long double>::infinity()) || CM::IsCloseToZero(diff, 0.001));
                // std::cout << std::fixed << std::setprecision(8) << res.results[m] << ", " << correctAreas[m] << std::endl;
                maxSolutionsDistance = std::max(maxSolutionsDistance, std::fabs(diff));
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
                EXPECT_TRUE((res.results[m].first == std::numeric_limits<long double>::infinity() &&
                             correctAreas[m] == std::numeric_limits<long double>::infinity()) || CM::IsCloseToZero(diff, 0.001));
                EXPECT_EQ(res.results[m].second, correctCounts[m]);
                maxSolutionsDistance = std::max(maxSolutionsDistance, std::fabs(diff));
            }
#ifdef VERBOSE
            const auto areThereCollinearPoints = AreThereCollinearPoints(points);
            std::cout   << std::fixed << std::setprecision(8) << "File " << (1+i) << "/" << (filesCount / 2) << ". Collinear points="
                        << (areThereCollinearPoints ? "yes" : "no") << ". Max distance=" << std::to_string(maxSolutionsDistance) << std::endl;
#endif
        }
    }
}

TEST(CrossCheckEppsteinAntipodal, BasicAssertions)
{
    constexpr auto shouldReconstructHull = true;
    const std::vector<std::string> paths = {
            "../../data/samples/eppstein/50",
            "../../data/samples/eppstein/60",
            "../../data/samples/eppstein/70",
            "../../data/samples/eppstein/80",
            "../../data/samples/eppstein/90",
            "../../data/samples/eppstein/100"
    };

    std::vector<long double> maximumAreas = {5, 10, 25, 50, 150, 200, 250, 300, 350, 400, std::numeric_limits<long double>::infinity()};
    std::vector<size_t> maximumCounts = {5, 10, 15, 20, 25, 30, 35, 40, 45, 50, (size_t) - 1};

    for (const auto & path : paths)
    {
        const auto filesCount = CountFilesInDirectory(path);
#ifdef VERBOSE
        std::cout << "TESTING SAMPLES IN FOLDER: " << path << std::endl;
#endif
        for(size_t i = 0; i < filesCount / 2; ++i)
        {
            const auto sampleFilename = std::string("data_") + std::to_string(i) + ".txt";
            std::vector<CM::Point2> points;
            SZ::ReadPointsFromFile(fs::path{path} / fs::path{sampleFilename}, points);

            for (const auto maxAllowedArea : maximumAreas)
            {
                for (const auto maxAllowedPointsCount : maximumCounts)
                {
                    const auto eppsteinRes = MT::EppsteinAlgorithm(points, maxAllowedPointsCount, maxAllowedArea, shouldReconstructHull);

                    std::vector<CM::Point2> eppsteinHull;
                    GetHullPoints(eppsteinRes.myHullIndices, points, eppsteinHull);
                    const auto eppsteinDiameter = ComputeDiameter(eppsteinHull);

                    constexpr auto epsilon = 0.00001l;
                    const auto maxAllowedDiameter = std::sqrtl(eppsteinDiameter.Value2()) + epsilon;

                    const auto antipodalRes = MT::AntipodalAlgorithm(points, maxAllowedPointsCount, maxAllowedArea, maxAllowedDiameter, shouldReconstructHull);

                    EXPECT_EQ(eppsteinRes.myHasFoundSolution, antipodalRes.myHasFoundSolution);
                    EXPECT_TRUE(antipodalRes.myHullIndices.size() == eppsteinRes.myHullIndices.size());

                    EXPECT_TRUE((eppsteinRes.myHullArea == std::numeric_limits<long double>::infinity() &&
                                 antipodalRes.myHullArea == std::numeric_limits<long double>::infinity()) ||
                                    CM::IsCloseToZero(eppsteinRes.myHullArea - antipodalRes.myHullArea, 0.001));

                    EXPECT_EQ(eppsteinRes.myPointsCount, antipodalRes.myPointsCount);
                    EXPECT_EQ(eppsteinDiameter, antipodalRes.myDiameter);

#ifdef VERBOSE
                    PrintResults(eppsteinRes, eppsteinDiameter, "Epp.Res");
                    PrintResults(antipodalRes, antipodalRes.myDiameter, "\tAnti.Res", true);
#endif
                }
            }



#ifdef VERBOSE
            const auto areThereCollinearPoints = AreThereCollinearPoints(points);
            std::cout   << std::fixed << std::setprecision(8) << "File " << (1+i) << "/" << (filesCount / 2) << ". Collinear points="
                        << (areThereCollinearPoints ? "yes" : "no") << std::endl;
#endif
        }
    }
}


// Cross check Eppstein and Antipodal from k=3 to n. First run Eppstein and store the area and diameter. Then use the area and diameter as constraints in Antipodal.
// Area and Diameter should match
