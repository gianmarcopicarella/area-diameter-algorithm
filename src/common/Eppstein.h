//
// Created by Gianmarco Picarella on 27/12/24.
//

#ifndef MASTER_THESIS_EPPSTEIN_H
#define MASTER_THESIS_EPPSTEIN_H

#include <vector>
#include "CustomMath.h"

#define DEBUG_EPPSTEIN

namespace MT
{
    struct EppsteinResult
    {
        bool myHasFoundSolution = false;
        float myHullArea = 0;
        int myPointsCount = 0;
        std::vector<int> myHullIndices;
#ifdef DEBUG_EPPSTEIN
        std::vector<long double> results;
#endif
    };

    EppsteinResult EppsteinAlgorithm(
            const std::vector<CM::Point2>& somePoints,
            int aMaxPointsCount,
            float aMaxArea = std::numeric_limits<float>::max(),
            bool aShouldReconstructHull = false);
}

#endif //MASTER_THESIS_EPPSTEIN_H
