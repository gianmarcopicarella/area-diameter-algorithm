#include "raylib.h"

#include "../common/Eppstein.h"
#include "../common/Parser.h"

#include <cassert>
#include <iostream>
using namespace MT;

int main(void)
{
    for(int i = 0; i < 30; ++i)
    {
        std::cout << i << std::endl;
        std::vector<CM::Point2> points;
        SZ::ReadPointsFromFile(std::string("../../data/samples/data_") + std::to_string(i) + ".txt", points);
        const auto res = MT::EppsteinAlgorithm(points, points.size());
        std::vector<long double> solutionsPY;
        SZ::ReadSolutionsFromFile(std::string("../../data/samples/results_") + std::to_string(i) + ".txt", solutionsPY);
        for(int m = 0; m < std::min(solutionsPY.size(), res.results.size()); ++m)
        {
            const auto temp = res.results[m] - solutionsPY[m];
            std::cout << temp << std::endl;
            if(std::isinf(temp) || !CM::IsCloseToZero(temp, 0.0001))
            {
                std::cout << "o" << std::endl;
            }
            assert(CM::IsCloseToZero(temp, 0.001));
        }
    }

    return 0;
}