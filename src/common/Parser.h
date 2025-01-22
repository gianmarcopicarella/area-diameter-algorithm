//
// Created by Gianmarco Picarella on 28/12/24.
//

#ifndef MASTER_THESIS_PARSER_H
#define MASTER_THESIS_PARSER_H

#include "CustomMath.h"

#include <string>
#include <fstream>
#include <cassert>
#include <vector>

namespace MT
{
    namespace SZ
    {
        bool ReadPointsFromFile(const std::string& aFilePath, std::vector<CM::Point2>& someOutPoints);

        template <typename T>
        bool ReadSolutionsFromFile(const std::string& aFilePath, std::vector<T>& someOutSolutions)
        {
            std::ifstream file(aFilePath);
            if(file.is_open())
            {
                int solutionsCount = 0;
                file >> solutionsCount;
                someOutSolutions.reserve(solutionsCount);
                T entry;
                while(file >> entry)
                {
                    someOutSolutions.emplace_back(entry);
                }
                assert(solutionsCount == someOutSolutions.size());
                return true;
            }
            return false;
        }
    }
}

#endif //MASTER_THESIS_PARSER_H
