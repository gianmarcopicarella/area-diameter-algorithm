//
// Created by Gianmarco Picarella on 06/01/25.
//

#ifndef MASTER_THESIS_ANTIPODAL_H
#define MASTER_THESIS_ANTIPODAL_H

#include <vector>
#include "CustomMath.h"

#define DEBUG_ANTIPODAL
// #define VERBOSE_ANTIPODAL

namespace MT
{
    struct AntipodalResult
    {
        long double myHullArea { std::numeric_limits<long double>::infinity() };
        size_t myPointsCount { 0 };
        bool myHasFoundSolution { false };
        std::vector<size_t> myHullIndices {};
#ifdef DEBUG_ANTIPODAL
        std::vector<std::pair<long double, size_t>> results;
#endif
    };

    AntipodalResult AntipodalAlgorithm(
            const std::vector<CM::Point2>& somePoints,
            size_t aMaxPointsCount,
            long double aMaxArea = std::numeric_limits<long double>::max(),
            long double aMaxDiameter = std::numeric_limits<long double>::max(),
            bool aShouldReconstructHull = false);
}


#endif //MASTER_THESIS_ANTIPODAL_H
