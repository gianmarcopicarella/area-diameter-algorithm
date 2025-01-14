//
// Created by Gianmarco Picarella on 27/12/24.
//

#ifndef MASTER_THESIS_EPPSTEIN_H
#define MASTER_THESIS_EPPSTEIN_H

#include <vector>
#include "CustomMath.h"

//#define DEBUG_EPPSTEIN

namespace MT
{
    struct EppsteinResult
    {
        bool myHasFoundSolution = false;
        long double myHullArea = 0;
        size_t myPointsCount = 0;
        std::vector<size_t> myHullIndices;
#ifdef DEBUG_EPPSTEIN
        std::vector<long double> results;
#endif
    };

    EppsteinResult EppsteinAlgorithm(
            const std::vector<CM::Point2>& somePoints,
            size_t aMaxPointsCount,
            long double aMaxArea = std::numeric_limits<long double>::max(),
            bool aShouldReconstructHull = false);
}

#endif //MASTER_THESIS_EPPSTEIN_H
