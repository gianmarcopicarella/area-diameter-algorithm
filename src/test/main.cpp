#include <gtest/gtest.h>

#include "../common/Eppstein.h"
#include "../common/Antipodal.h"
#include "../common/Parser.h"

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

TEST(EppsteinWithSolutions, BasicAssertions)
{
    const std::vector<std::string> paths = {
            "../../data/samples/eppstein/10",
            "../../data/samples/eppstein/20",
            "../../data/samples/eppstein/30",
            "../../data/samples/eppstein/40",
            "../../data/samples/eppstein/50",
            "../../data/samples/eppstein/60",
            "../../data/samples/eppstein/70",
            "../../data/samples/eppstein/80",
            "../../data/samples/eppstein/90",
            "../../data/samples/eppstein/100"
    };

    for (const auto & path : paths)
    {
#ifdef VERBOSE
        std::cout << "TESTING SAMPLES IN FOLDER: " << path << std::endl;
#endif
        std::vector<SZ::Solution> solutions;
        SZ::ReadSolutionsFromFile(fs::path{path} / fs::path{"results.csv"}, solutions);
        for (size_t i = 0; i < CountPointFilesInFolder(path); ++i)
        {
            std::vector<CM::Point2> points;
            const auto filename = std::string("points_") + std::to_string(i) + ".csv";
            SZ::ReadPointsFromFile(fs::path{path} / fs::path{filename}, points);

            const auto maxArea = std::get<0>(solutions[i]);
            const auto maxPointsCount = std::get<2>(solutions[i]);
            const auto convexArea = std::get<3>(solutions[i]);

            constexpr auto reconstructHull = true;
            const auto result = MT::EppsteinAlgorithm(points, maxPointsCount, maxArea, reconstructHull);

            EXPECT_EQ(result.has_value(), convexArea.has_value());
            if (result.has_value())
            {
                EXPECT_EQ(result->myPointsCount, convexArea->myPointsCount);
                EXPECT_EQ(result->myDiameterOpt, convexArea->myDiameterOpt);
                EXPECT_FLOAT_EQ(result->myHullArea, convexArea->myHullArea);
            }

#ifdef VERBOSE
            const auto areThereCollinearPoints = AreThereCollinearPoints(points);
            std::cout << std::fixed << std::setprecision(8) << "File " << filename << ". Collinear points="
                      << (areThereCollinearPoints ? "yes" : "no") << std::endl;
#endif
        }
    }
}

TEST(AntipodalWithEppstein, BasicAssertions)
{
    const std::vector<std::string> paths = {
            "../../data/samples/eppstein/10",
            "../../data/samples/eppstein/20",
            "../../data/samples/eppstein/30",
            "../../data/samples/eppstein/40",
            "../../data/samples/eppstein/50",
            "../../data/samples/eppstein/60",
            "../../data/samples/eppstein/70",
            "../../data/samples/eppstein/80",
            "../../data/samples/eppstein/90",
            "../../data/samples/eppstein/100"
    };

    constexpr auto epsilon = 0.000001l;
    constexpr auto reconstructHull = true;

    for (const auto & path : paths)
    {
#ifdef VERBOSE
        std::cout << "TESTING SAMPLES IN FOLDER: " << path << std::endl;
#endif
        std::vector<SZ::Solution> solutions;
        SZ::ReadSolutionsFromFile(fs::path{path} / fs::path{"results.csv"}, solutions);
        for (size_t i = 0; i < CountPointFilesInFolder(path); ++i)
        {
            std::vector<CM::Point2> points;
            const auto filename = std::string("points_") + std::to_string(i) + ".csv";
            SZ::ReadPointsFromFile(fs::path{path} / fs::path{filename}, points);

            const auto maxArea = std::get<0>(solutions[i]);
            const auto maxPointsCount = std::get<2>(solutions[i]);
            const auto convexArea = std::get<3>(solutions[i]);

            const auto result = MT::EppsteinAlgorithm(points, maxPointsCount, maxArea, reconstructHull);

            EXPECT_EQ(result.has_value(), convexArea.has_value());
            if (result.has_value())
            {
                EXPECT_EQ(result->myPointsCount, convexArea->myPointsCount);
                EXPECT_EQ(result->myDiameterOpt, convexArea->myDiameterOpt);
                EXPECT_FLOAT_EQ(result->myHullArea, convexArea->myHullArea);

                const auto maxDiameter = std::sqrtl(
                        CM::Distance2(points[result->myDiameterOpt->myFirstIndex], points[result->myDiameterOpt->mySecondIndex]));
                const auto antipodalResult = MT::AntipodalAlgorithm(points, maxPointsCount, maxArea, maxDiameter + epsilon, reconstructHull);

                EXPECT_EQ(result->myPointsCount, antipodalResult->myPointsCount);
                EXPECT_EQ(result->myDiameterOpt, antipodalResult->myDiameterOpt);
                EXPECT_FLOAT_EQ(result->myHullArea, antipodalResult->myHullArea);
            }
            else
            {
                constexpr auto maxDiameter = std::numeric_limits<long double>::infinity();
                const auto antipodalResult = MT::AntipodalAlgorithm(points, maxPointsCount, maxArea, maxDiameter, reconstructHull);
                EXPECT_EQ(result.has_value(), antipodalResult.has_value());
            }

#ifdef VERBOSE
            const auto areThereCollinearPoints = AreThereCollinearPoints(points);
            std::cout << std::fixed << std::setprecision(8) << "File " << filename << ". Collinear points="
                      << (areThereCollinearPoints ? "yes" : "no") << std::endl;
#endif
        }
    }
}
