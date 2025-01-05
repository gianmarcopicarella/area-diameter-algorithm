//
// Created by Gianmarco Picarella on 27/12/24.
//

#ifndef MASTER_THESIS_CUSTOMMATH_H
#define MASTER_THESIS_CUSTOMMATH_H

#include <limits>
#include <vector>

namespace MT
{
    namespace CM
    {
        struct Point2
        {
            Point2() = default;
            Point2(long double x, long double y, size_t index) :
                x(x), y(y), index(index) {}
            long double x, y;
            size_t index;
        };

        constexpr auto SortPointsHorizontally = [](const Point2& p1, const Point2& p2){
            return (p1.x < p2.x) || (p1.x == p2.x && p1.y < p2.y);
        };

        constexpr auto SortPointsVertically = [](const Point2& p1, const Point2& p2){
            return p1.y < p2.y;
        };

        enum class ORIENTATION
        {
            RIGHT_TURN = -1,
            COLLINEAR = 0,
            LEFT_TURN = 1
        };

        long double Dot2(const Point2& aFirstPoint, const Point2& aSecondPoint);

        long double SquaredDistance(const Point2& aFirstPoint, const Point2& aSecondPoint);

        bool IsCloseToZero(long double aValue,
                           long double anEpsilon = std::numeric_limits<long double>::epsilon());

        bool AreCollinear(const Point2& aFirstPoint, const Point2& aSecondPoint, const Point2& aThirdPoint);

        long double SignedArea(const Point2& aFirstPoint, const Point2& aSecondPoint, const Point2& aThirdPoint);

        ORIENTATION Orientation(const Point2& aFirstPoint, const Point2& aSecondPoint, const Point2& aThirdPoint);

        long double Angle(const CM::Point2& aReferencePoint, const CM::Point2& aPoint);
    }
}

#endif //MASTER_THESIS_CUSTOMMATH_H
