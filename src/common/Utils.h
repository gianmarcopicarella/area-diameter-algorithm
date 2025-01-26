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
    struct Diameter
    {
        CM::Point2 myFirstPoint {}, mySecondPoint {};
        long double Value2() const;
        bool operator==(const Diameter&) const;
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

    Diameter ComputeDiameter(const std::vector<CM::Point2>& somePoints);

    template<typename Container>
    void ForAllAntipodalPairs(const Container& somePoints,
                              const std::function<void(const CM::Point2&, const CM::Point2&)>& aFunction,
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
            aFunction(somePoints[pi], somePoints[qi]);

            while(  CM::SignedArea(somePoints[pi], somePoints[NEXT(pi, pointsCount)], somePoints[NEXT(qi, pointsCount)]) >
                    CM::SignedArea(somePoints[pi], somePoints[NEXT(pi, pointsCount)], somePoints[qi]))
            {
                qi = NEXT(qi, pointsCount);
                if(pi != q0 || qi != 0)
                {
                    aFunction(somePoints[pi], somePoints[qi]);
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
                    aFunction(somePoints[pi], somePoints[NEXT(qi, pointsCount)]);
                }
                else
                {
                    break;
                }
            }
        }
    }
}


#endif //MASTER_THESIS_UTILS_H
