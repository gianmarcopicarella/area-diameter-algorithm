//
// Created by Gianmarco Picarella on 28/12/24.
//

#include "Parser.h"

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
    }
}