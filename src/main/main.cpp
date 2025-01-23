#include "../common/Antipodal.h"
#include "../common/Parser.h"
#include "../common/Utils.h"

#include <cassert>
#include <iostream>
#include <iomanip>
#include <chrono>

#include <filesystem>
namespace fs = std::filesystem;

using namespace MT;

int main(void)
{
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::seconds;

    std::vector<CM::Point2> points;
    SZ::ReadPointsFromFile(fs::path{"../../data/samples/antipodal/100/data_0.txt"}, points);

    constexpr auto maxDiameter = std::numeric_limits<long double>::infinity();
    constexpr auto maxArea = std::numeric_limits<long double>::infinity();
    constexpr auto reconstructHull = true;

    auto t1 = high_resolution_clock::now();
    const auto res = AntipodalAlgorithm(points, points.size(), maxArea, maxDiameter, reconstructHull);
    auto t2 = high_resolution_clock::now();

    auto sec_int = duration_cast<seconds>(t2 - t1);
    std::cout << "The program took " << sec_int.count() << " seconds to run" << std::endl;
    for(auto i : res.myHullIndices) std::cout << i << ", ";
    std::cout << std::endl;

    return 0;
}