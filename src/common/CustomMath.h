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
            Point2(long double anX, long double anY, size_t anIndex) :
                    myX(anX), myY(anY), myIndex(anIndex) {}
            long double myX, myY;
            size_t myIndex;
            bool operator==(const Point2&) const;
        };

        constexpr auto SortPointsHorizontally = [](const Point2& p1, const Point2& p2){
            return (p1.myX < p2.myX) || (p1.myX == p2.myX && p1.myY < p2.myY);
        };

        constexpr auto SortPointsVertically = [](const Point2& p1, const Point2& p2){
            return p1.myY < p2.myY;
        };

        enum class ORIENTATION
        {
            RIGHT_TURN = -1,
            COLLINEAR = 0,
            LEFT_TURN = 1
        };

        template <typename T> int Sign(const T& aValue)
        {
            return (T(0) < aValue) - (aValue < T(0));
        }

        long double Dot(const Point2& aFirstPoint, const Point2& aSecondPoint);

        long double Distance2(const Point2& aFirstPoint, const Point2& aSecondPoint);

        // Using this value for epsilon guarantees that std::atan2l and SignedArea produce coherent results for collinear points
        bool IsCloseToZero(long double aValue,
                           long double anEpsilon = 0.0000000001l/*std::numeric_limits<long double>::epsilon()*/);

        bool AreCollinear(const Point2& aFirstPoint, const Point2& aSecondPoint, const Point2& aThirdPoint);

        long double SignedArea(const Point2& aFirstPoint, const Point2& aSecondPoint, const Point2& aThirdPoint);

        ORIENTATION Orientation(const Point2& aFirstPoint, const Point2& aSecondPoint, const Point2& aThirdPoint);

        long double Angle(const CM::Point2& aReferencePoint, const CM::Point2& aPoint);

        long double ProjectedDistance2(const CM::Point2& aStartPoint, const CM::Point2& anEndPoint, const CM::Point2& aPoint);

        long double PointLinePseudoDistance(const CM::Point2& aStartPoint, const CM::Point2& anEndPoint, const CM::Point2& aPoint);
    }
}

#endif //MASTER_THESIS_CUSTOMMATH_H
