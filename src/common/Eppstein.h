#ifndef MASTER_THESIS_EPPSTEIN_H
#define MASTER_THESIS_EPPSTEIN_H

#include <optional>

#include "CustomMath.h"
#include "Utils.h"

namespace MT
{
    std::optional<ConvexArea> EppsteinAlgorithm(
            const std::vector<CM::Point2>& somePoints,
            size_t aMaxPointsCount = (size_t) - 1,
            long double aMaxArea = std::numeric_limits<long double>::max(),
            bool aShouldReconstructHull = false,
            bool aShouldEnableOptimizations = false);

    std::optional<ConvexArea> EppsteinAlgorithmWithBenchmarkInfo(
            const std::vector<CM::Point2>& somePoints,
            std::optional<BenchmarkInfo>& anOutBenchmarkInfoOpt,
            size_t aMaxPointsCount = (size_t) - 1,
            long double aMaxArea = std::numeric_limits<long double>::max(),
            bool aShouldReconstructHull = false,
            bool aShouldEnableOptimizations = false);
}

#endif //MASTER_THESIS_EPPSTEIN_H
