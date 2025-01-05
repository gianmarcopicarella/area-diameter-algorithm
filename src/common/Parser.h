//
// Created by Gianmarco Picarella on 28/12/24.
//

#ifndef MASTER_THESIS_PARSER_H
#define MASTER_THESIS_PARSER_H

#include "CustomMath.h"

#include <string>

namespace MT
{
    namespace SZ
    {
        bool ReadPointsFromFile(const std::string& aFilePath, std::vector<CM::Point2>& someOutPoints);
        bool ReadSolutionsFromFile(const std::string& aFilePath, std::vector<long double>& someOutSolutions);
    }
}

#endif //MASTER_THESIS_PARSER_H
