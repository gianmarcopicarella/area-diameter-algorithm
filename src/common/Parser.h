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
#include <optional>

namespace MT
{
    namespace SZ
    {
        bool ReadPointsFromFile(const std::string& aFilePath, std::vector<CM::Point2>& someOutPoints);

        using Metadata = std::optional<std::tuple<long double, long double, size_t>>;
        template<typename... U>
        Metadata ReadSolutionsFromFile(const std::string& aFilePath, std::vector<U>& ... someOutArgs)
        {
            static_assert(sizeof...(U) == sizeof...(someOutArgs));
            std::ifstream file(aFilePath);
            Metadata metadata;
            if(file.is_open())
            {
                size_t entriesCount;
                file >> entriesCount;
                entriesCount /= sizeof...(someOutArgs);
                {
                    long double maxArea, maxDiameter;
                    size_t maxCount;
                    file >> maxArea;
                    file >> maxDiameter;
                    file >> maxCount;
                    metadata = {maxArea, maxDiameter, maxCount};
                }
                ([&](auto& vec) {
                    using ValueType = typename std::decay_t<decltype(vec)>::value_type;
                    ValueType value;
                    for(size_t i = 0; i < entriesCount; ++i)
                    {
                        file >> value;
                        vec.push_back(value);
                    }
                }(someOutArgs), ...);
            }
            return metadata;
        }
    }
}

#endif //MASTER_THESIS_PARSER_H
