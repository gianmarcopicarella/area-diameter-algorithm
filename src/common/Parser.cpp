//
// Created by Gianmarco Picarella on 28/12/24.
//

#include "Parser.h"
#include <fstream>
#include <assert.h>
#include <vector>

namespace MT
{
    namespace SZ
    {
        bool ReadPointsFromFile(const std::string& aFilePath, std::vector<CM::Point2>& someOutPoints)
        {
            std::ifstream file(aFilePath);
            if(file.is_open())
            {
                int pointsCount = 0;
                file >> pointsCount;
                someOutPoints.reserve(pointsCount);
                long double x, y;
                size_t index = 0;
                while(file >> x >> y)
                {
                    someOutPoints.emplace_back(x, y, index++);
                }
                assert(pointsCount == someOutPoints.size());
                return true;
            }
            return false;
        }

        bool ReadSolutionsFromFile(const std::string& aFilePath, std::vector<long double>& someOutSolutions)
        {
            std::ifstream file(aFilePath);
            if(file.is_open())
            {
                int solutionsCount = 0;
                file >> solutionsCount;
                someOutSolutions.reserve(solutionsCount);
                long double minArea;
                while(file >> minArea)
                {
                    someOutSolutions.emplace_back(minArea);
                }
                assert(solutionsCount == someOutSolutions.size());
                return true;
            }
            return false;
        }
    }
}