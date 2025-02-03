//
// Created by Gianmarco Picarella on 06/01/25.
//

#ifndef MASTER_THESIS_ANTIPODAL_H
#define MASTER_THESIS_ANTIPODAL_H

#include <vector>
#include "CustomMath.h"
#include "Utils.h"

// #define OPT_PRE_SORT_SEGMENTS
#define OPT_USE_OPTIMAL_SOLUTION

namespace MT
{
    std::optional<ConvexArea> AntipodalAlgorithm(
            const std::vector<CM::Point2>& somePoints,
            size_t aMaxPointsCount = (size_t) - 1,
            long double aMaxArea = std::numeric_limits<long double>::max(),
            long double aMaxDiameter = std::numeric_limits<long double>::max(),
            bool aShouldReconstructHull = false);
}


#endif //MASTER_THESIS_ANTIPODAL_H
