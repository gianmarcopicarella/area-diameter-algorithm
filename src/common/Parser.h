#ifndef MASTER_THESIS_PARSER_H
#define MASTER_THESIS_PARSER_H

#include "CustomMath.h"
#include "Utils.h"

#include <string>
#include <fstream>
#include <cassert>
#include <vector>
#include <optional>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

namespace MT
{
    struct Solution
    {
        std::string myName { "None" };
        long double myMaxArea { std::numeric_limits<long double>::infinity() };
        long double myMaxDiameter { std::numeric_limits<long double>::infinity() };
        size_t myMaxCount { (size_t) - 1 };
        std::optional<ConvexArea> myConvexAreaOpt;
    };

    void from_json(const json& aJson, Solution& anOutSolution);
    void to_json(json& aJson, const Solution& aSolution);

    namespace SZ
    {
        void ReadPointsFromFile(const std::string& aFilePath, std::vector<CM::Point2>& someOutPoints);
        void ReadSolutionsFromFile(const std::string& aFilePath, std::vector<Solution>& someOutSolutions);
    }
}

#endif //MASTER_THESIS_PARSER_H
