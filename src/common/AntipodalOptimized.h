//
// Created by Gianmarco Picarella on 30/03/25.
//

#ifndef MASTER_THESIS_ANTIPODALOPTIMIZED_H
#define MASTER_THESIS_ANTIPODALOPTIMIZED_H

#include <vector>
#include "CustomMath.h"
#include "Utils.h"

namespace MT
{
    std::optional<ConvexArea> AntipodalOptimizedAlgorithm(
            const std::vector<CM::Point2>& somePoints,
            size_t aMaxPointsCount = (size_t) - 1,
            long double aMaxArea = std::numeric_limits<long double>::max(),
            long double aMaxDiameter = std::numeric_limits<long double>::max(),
            bool aShouldReconstructHull = false,
            bool aShouldEnableOptimizations = false);

    std::optional<ConvexArea> AntipodalOptimizedAlgorithmWithBenchmarkInfo(
            const std::vector<CM::Point2>& somePoints,
            std::optional<BenchmarkInfo>& anOutBenchmarkInfoOpt,
            size_t aMaxPointsCount = (size_t) - 1,
            long double aMaxArea = std::numeric_limits<long double>::max(),
            long double aMaxDiameter = std::numeric_limits<long double>::max(),
            bool aShouldReconstructHull = false,
            bool aShouldEnableOptimizations = false);
}

#endif //MASTER_THESIS_ANTIPODALOPTIMIZED_H
