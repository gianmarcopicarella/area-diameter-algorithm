//
// Created by Gianmarco Picarella on 28/12/24.
//

#ifndef MASTER_THESIS_PARSER_H
#define MASTER_THESIS_PARSER_H

#include "CustomMath.h"
#include "Utils.h"

#include <string>
#include <fstream>
#include <cassert>
#include <vector>
#include <optional>

namespace MT
{
    namespace SZ
    {
        void ReadPointsFromFile(const std::string& aFilePath, std::vector<CM::Point2>& someOutPoints);
        using Solution = std::tuple<long double, long double, size_t, std::optional<ConvexArea>>;
        void ReadSolutionsFromFile(const std::string& aFilePath, std::vector<Solution>& someOutSolutions);
        void WriteLineToCsv(const std::string& aFilePath, const std::function<std::string(void)>& aCsvLineGenerator, bool aShouldAppendToFile = true);
    }
}

#endif //MASTER_THESIS_PARSER_H
