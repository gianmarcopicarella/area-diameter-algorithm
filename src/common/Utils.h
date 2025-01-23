//
// Created by Gianmarco Picarella on 28/12/24.
//

#ifndef MASTER_THESIS_UTILS_H
#define MASTER_THESIS_UTILS_H

#include "CustomMath.h"
#include <vector>

namespace MT
{
    void CountPointsBelowAllSegments(const std::vector<CM::Point2>& somePoints,
                                     const std::vector<std::vector<CM::Point2>>& someClockWiseSortedPoints,
                                     std::vector<std::vector<int>>& someOutBelowCounts,
                                     std::vector<std::vector<int>>& someOutCollinearCounts);

    int PointsInTriangle(
            const CM::Point2& aFirstPoint,
            const CM::Point2& aSecondPoint,
            const CM::Point2& aThirdPoint,
            const std::vector<std::vector<int>>& somePointCountBelowSegments,
            const std::vector<std::vector<int>>& someCollinearPointCounts);

    bool ArePointsClockwise(const CM::Point2& aReferencePoint, const CM::Point2& aFirstPoint, const CM::Point2& aSecondPoint);

    void SortPointsClockWiseAroundPoint(const CM::Point2& aReferencePoint, std::vector<CM::Point2>& someOutClockWiseSortedPoints);

    long double ComputeDiameter(const std::vector<CM::Point2>& somePoints);

    void FindAntipodalPairs(
            const std::vector<CM::Point2>& somePoints,
            std::vector<std::pair<size_t, size_t>>& someOutAntipodalIndices);

    void AreAntipodalPairs( const std::array<CM::Point2, 6>& somePoints,
                            size_t aLen,
                            size_t aFirstIndex, size_t aSecondIndex,
                            size_t aThirdIndex, size_t aFourthIndex,
                            bool& anOutIsFirstAntipodal, bool& anOutIsSecondAntipodal);
}


#endif //MASTER_THESIS_UTILS_H
