#include <gtest/gtest.h>

#include "../common/Eppstein.h"
#include "../common/Antipodal.h"
#include "../common/AntipodalOptimized.h"
#include "../common/Parser.h"
#include "../common/Constants.h"

#define VERBOSE

#include <iostream>
#include <filesystem>
namespace fs = std::filesystem;

using namespace MT;

size_t CountPointFilesInFolder(const std::string& aPath)
{
    size_t count = 0;
    for (const auto& entry : fs::directory_iterator(aPath))
    {
        if(entry.path().filename().string().find("points_") != std::string::npos)
        {
            ++count;
        }
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

void CheckHullIndices(const std::vector<size_t> & someIndices, const std::vector<size_t> & someOtherIndices)
{
    EXPECT_TRUE((someIndices.empty() && someOtherIndices.empty()) ||
                (someIndices.size() == someOtherIndices.size()));
    if(!someIndices.empty())
    {
        const auto referenceIndexIter =
                std::find(someIndices.begin(), someIndices.end(), someOtherIndices[0]);
        EXPECT_TRUE(referenceIndexIter != someIndices.end());
        const auto referenceIndex = std::distance(someIndices.begin(), referenceIndexIter);
        for(size_t i = 0, j = referenceIndex; i < someIndices.size(); ++i, j = ((j + 1) % someIndices.size()))
        {
            EXPECT_EQ(someOtherIndices[i], someIndices[j]);
        }
    }
}

TEST(CrossTestEppsteinAndAntipodal, BasicAssertions)
{
    for (const auto& entry : fs::directory_iterator(MT::Constants::PATH_TO_TESTS))
    {
        if(!entry.is_directory())
        {
            continue;
        }

#ifdef VERBOSE
        std::cout << "Running tests in folder [" << entry.path() << "]" << std::endl;
#endif
        std::vector<Solution> solutions;
        SZ::ReadSolutionsFromFile(entry.path() / fs::path{"results.json"}, solutions);
        for (size_t i = 0; i < CountPointFilesInFolder(entry.path()); ++i)
        {
            std::vector<CM::Point2> points;
            const auto filename = std::string("points_") + std::to_string(i) + ".json";
            SZ::ReadPointsFromFile(entry.path() / fs::path{filename}, points);

            const auto& convexArea = solutions[i].myConvexAreaOpt;

            constexpr auto reconstructHull = true;
            constexpr auto epsilon = 0.000001l;

            const auto result = MT::EppsteinAlgorithm(points, solutions[i].myMaxCount, solutions[i].myMaxArea, reconstructHull);

            EXPECT_EQ(result.has_value(), convexArea.has_value());
            if (result.has_value())
            {
                EXPECT_EQ(result->myPointsCount, convexArea->myPointsCount);
                EXPECT_EQ(result->myDiameterOpt, convexArea->myDiameterOpt);
                EXPECT_FLOAT_EQ(result->myHullArea, convexArea->myHullArea);
                CheckHullIndices(result->myHullIndices, convexArea->myHullIndices);

                const auto maxDiameter = std::sqrtl(
                        CM::Distance2(points[result->myDiameterOpt->myFirstIndex], points[result->myDiameterOpt->mySecondIndex]));
                const auto antipodalResult = MT::AntipodalAlgorithm(points, solutions[i].myMaxCount, solutions[i].myMaxArea, maxDiameter + epsilon, reconstructHull);

                EXPECT_EQ(result->myPointsCount, antipodalResult->myPointsCount);
                EXPECT_EQ(result->myDiameterOpt, antipodalResult->myDiameterOpt);
                EXPECT_FLOAT_EQ(result->myHullArea, antipodalResult->myHullArea);
                CheckHullIndices(result->myHullIndices, antipodalResult->myHullIndices);

                // Optimized antipodal algorithm
                const auto optAntipodalResult = MT::AntipodalOptimizedAlgorithm(points, solutions[i].myMaxCount, solutions[i].myMaxArea, maxDiameter + epsilon, reconstructHull);
                EXPECT_EQ(result->myPointsCount, optAntipodalResult->myPointsCount);
                EXPECT_EQ(result->myDiameterOpt, optAntipodalResult->myDiameterOpt);
                EXPECT_FLOAT_EQ(result->myHullArea, optAntipodalResult->myHullArea);
                // Missing hull indices test
            }
            else
            {
                constexpr auto maxDiameter = std::numeric_limits<long double>::infinity();
                const auto antipodalResult = MT::AntipodalAlgorithm(points, solutions[i].myMaxCount, solutions[i].myMaxArea, maxDiameter, reconstructHull);
                EXPECT_EQ(result.has_value(), antipodalResult.has_value());

                // Optimized antipodal algorithm
                // const auto optAntipodalResult = MT::AntipodalOptimizedAlgorithm(points, solutions[i].myMaxCount, solutions[i].myMaxArea, maxDiameter, reconstructHull);
                // EXPECT_EQ(result.has_value(), optAntipodalResult.has_value());
            }

#ifdef VERBOSE
            const auto areThereCollinearPoints = AreThereCollinearPoints(points);
            std::cout << "Done with file [" << filename << "]; Collinear points ["
                      << (areThereCollinearPoints ? "yes]" : "no]") << std::endl;
#endif
        }
    }
}
