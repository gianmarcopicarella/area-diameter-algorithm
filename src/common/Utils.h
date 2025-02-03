//
// Created by Gianmarco Picarella on 28/12/24.
//

#ifndef MASTER_THESIS_UTILS_H
#define MASTER_THESIS_UTILS_H

#include "CustomMath.h"
#include <vector>
#include <functional>
#include <cassert>

#define NEXT(x, n) ((x + 1) % (n))

namespace MT
{
    constexpr auto INVALID_INDEX = (size_t) - 1;
    struct Diameter
    {
        bool operator==(const Diameter&) const;
        size_t myFirstIndex { INVALID_INDEX }, mySecondIndex { INVALID_INDEX };
    };

    struct ConvexArea
    {
        long double myHullArea { std::numeric_limits<long double>::infinity() };
        size_t myPointsCount { 0 };
        std::optional<Diameter> myDiameterOpt;
        std::vector<size_t> myHullIndices {};
    };

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

    std::optional<Diameter> ComputeDiameter(const std::vector<CM::Point2>& somePoints);

    template<typename Container>
    void ForAllAntipodalPairs(const Container& somePoints,
                              const std::function<void(const size_t, const size_t)>& aFunction,
                              size_t aMaxLength = (size_t) - 1)
    {
        assert(somePoints.size() > 2);
        const auto pointsCount = aMaxLength == ((size_t) - 1) ? somePoints.size() : aMaxLength;
        size_t pi = pointsCount - 1;
        size_t qi = 0;

        while(  CM::SignedArea(somePoints[pi], somePoints[NEXT(pi, pointsCount)], somePoints[NEXT(qi, pointsCount)]) >
                CM::SignedArea(somePoints[pi], somePoints[NEXT(pi, pointsCount)], somePoints[qi]))
        {
            qi = NEXT(qi, pointsCount);
        }

        auto q0 = qi;
        while (qi != 0)
        {
            pi = NEXT(pi, pointsCount);
            aFunction(pi, qi);

            while(  CM::SignedArea(somePoints[pi], somePoints[NEXT(pi, pointsCount)], somePoints[NEXT(qi, pointsCount)]) >
                    CM::SignedArea(somePoints[pi], somePoints[NEXT(pi, pointsCount)], somePoints[qi]))
            {
                qi = NEXT(qi, pointsCount);
                if(pi != q0 || qi != 0)
                {
                    aFunction(pi, qi);
                }
                else
                {
                    break;
                }
            }

            const auto signedAreaDifference =
                    CM::SignedArea(somePoints[pi], somePoints[NEXT(pi, pointsCount)], somePoints[NEXT(qi, pointsCount)]) -
                    CM::SignedArea(somePoints[pi], somePoints[NEXT(pi, pointsCount)], somePoints[qi]);

            if(CM::IsCloseToZero(signedAreaDifference))
            {
                if(pi != q0 || qi != (pointsCount - 1))
                {
                    aFunction(pi, NEXT(qi, pointsCount));
                }
                else
                {
                    break;
                }
            }
        }
    }

    void GetHullPoints(const std::vector<size_t>& someIndices,
                       const std::vector<CM::Point2>& somePoints,
                       std::vector<CM::Point2>& someOutHullPoints);
}


#endif //MASTER_THESIS_UTILS_H
