#include "../common/Antipodal.h"
#include "../common/Eppstein.h"
#include "../common/Parser.h"
#include "../common/Utils.h"

#include <cassert>
#include <iostream>
#include <iomanip>
#include <chrono>

#include <filesystem>
namespace fs = std::filesystem;

using namespace MT;

void GetHullPoints(const std::vector<size_t>& someIndices,
                   const std::vector<CM::Point2>& somePoints,
                   std::vector<CM::Point2>& someOutHullPoints)
{
    someOutHullPoints.resize(someIndices.size());
    std::transform(someIndices.begin(), someIndices.end(), someOutHullPoints.begin(),
                   [&](auto index){ return somePoints[index]; });
}

int main(void)
{
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::seconds;

    for(size_t i = 1; i <= 10; ++i)
    {
        std::vector<CM::Point2> points;
        SZ::ReadPointsFromFile(fs::path{"../../data/samples/experiments/data_std_" + std::to_string(i) + ".txt"}, points);
        constexpr auto maxDiameter = std::numeric_limits<long double>::infinity();
        constexpr auto maxArea = 4.l;
        constexpr auto reconstructHull = true;

        auto t1 = high_resolution_clock::now();
        const auto res = AntipodalAlgorithm(points, points.size(), maxArea, maxDiameter, reconstructHull);
        auto t2 = high_resolution_clock::now();

        auto t3 = high_resolution_clock::now();
        const auto res2 = EppsteinAlgorithm(points, points.size(), maxArea, reconstructHull);
        auto t4 = high_resolution_clock::now();

        std::cout << "[" << res.myHullArea << ", " << res.myPointsCount << ", " << std::sqrtl(res.myDiameter.Value2()) << ", " << duration_cast<std::chrono::milliseconds>(t2 - t1).count() << ", [";
        for(int x = 0; x < res.myHullIndices.size(); ++x)
        {
            if(x < res.myHullIndices.size() - 1)
            {
                std::cout << res.myHullIndices[x] << ", ";
            }
            else
            {
                std::cout << res.myHullIndices[x];
            }
        }
        std::cout << "]]" << std::endl;

        std::vector<CM::Point2> eppHull;
        GetHullPoints(res2.myHullIndices, points, eppHull);

        std::cout << "[" << res2.myHullArea << ", " << res2.myPointsCount << ", " << std::sqrtl(ComputeDiameter(eppHull).Value2()) << ", " << duration_cast<std::chrono::milliseconds>(t4 - t3).count() << ", [";
        for(int x = 0; x < res2.myHullIndices.size(); ++x)
        {
            if(x < res2.myHullIndices.size() - 1)
            {
                std::cout << res2.myHullIndices[x] << ", ";
            }
            else
            {
                std::cout << res2.myHullIndices[x];
            }
        }
        std::cout << "]]" << std::endl;


/*
        std::cout << "Antipodal data:\n";
        std::cout << "area:" << res.myHullArea << ", count:" << res.myPointsCount << ", diam:" << std::sqrtl(res.myDiameter.Value2()) << std::endl;
        std::cout << "Time (ms):" << duration_cast<std::chrono::milliseconds>(t2 - t1).count() << std::endl;
        for(auto idx : res.myHullIndices) std::cout << idx << ", ";
        std::cout << std::endl << std::endl;

        std::cout << "Eppstein data:\n";
        std::cout << "area:" << res2.myHullArea << ", count:" << res2.myPointsCount << ", diam:" << std::sqrtl(ComputeDiameter(eppHull).Value2()) << std::endl;
        std::cout << "Time (ms):" << duration_cast<std::chrono::milliseconds>(t4 - t3).count() << std::endl;
        for(auto idx : res2.myHullIndices) std::cout << idx << ", ";
        std::cout << std::endl << std::endl;
        */
    }

    return 0;
}