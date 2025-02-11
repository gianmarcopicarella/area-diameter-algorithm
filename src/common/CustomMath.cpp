//
// Created by Gianmarco Picarella on 27/12/24.
//

#include "CustomMath.h"

#include <cmath>
#include <cassert>

namespace MT
{
    namespace CM
    {
        bool Point2::operator==(const Point2& anotherPoint) const
        {
            return anotherPoint.myX == myX && anotherPoint.myY == myY && anotherPoint.myIndex == myIndex;
        }

        long double Distance2(const Point2& aFirstPoint, const Point2& aSecondPoint)
        {
            const long double dx = aSecondPoint.myX - aFirstPoint.myX;
            const long double dy = aSecondPoint.myY - aFirstPoint.myY;
            return dx * dx + dy * dy;
        }

        bool IsCloseToZero(long double aValue, long double anEpsilon)
        {
            return std::fabsl(aValue) < anEpsilon;
        }

        long double Dot(const Point2& aFirstPoint, const Point2& aSecondPoint)
        {
            return aFirstPoint.myX * aSecondPoint.myX + aFirstPoint.myY * aSecondPoint.myY;
        }

        long double SignedArea(const Point2& aFirstPoint, const Point2& aSecondPoint, const Point2& aThirdPoint)
        {
            return 0.5l * ( (aSecondPoint.myX - aFirstPoint.myX) * (aThirdPoint.myY - aFirstPoint.myY) -
                            (aThirdPoint.myX - aFirstPoint.myX) * (aSecondPoint.myY - aFirstPoint.myY));
        }

        bool AreCollinear(const Point2& aFirstPoint, const Point2& aSecondPoint, const Point2& aThirdPoint)
        {
            return IsCloseToZero(SignedArea(aFirstPoint, aSecondPoint, aThirdPoint));
        }

        ORIENTATION Orientation(const Point2& aFirstPoint, const Point2& aSecondPoint, const Point2& aThirdPoint)
        {
            const auto signedArea = SignedArea(aFirstPoint, aSecondPoint, aThirdPoint);
            if(signedArea > 0)
            {
                return ORIENTATION::LEFT_TURN;
            }
            else if(signedArea < 0)
            {
                return ORIENTATION::RIGHT_TURN;
            }
            else
            {
                // assert(IsCloseToZero(signedArea));
                return ORIENTATION::COLLINEAR;
            }
        }

        long double Angle(const CM::Point2& aReferencePoint, const CM::Point2& aPoint)
        {
            const auto angle =
                    std::atan2l(aPoint.myY - aReferencePoint.myY, aPoint.myX - aReferencePoint.myX);
            if(angle < 0.f)
            {
                return angle + 2.f * M_PI;
            }
            return angle;
        }

        long double ProjectedDistance2(const CM::Point2& aStartPoint, const CM::Point2& anEndPoint, const CM::Point2& aPoint)
        {
            if(aPoint.myIndex == aStartPoint.myIndex)
            {
                return 0;
            }
            const auto dist2 = CM::Distance2(aStartPoint, anEndPoint);
            if(aPoint.myIndex == anEndPoint.myIndex)
            {
                return dist2;
            }
            else
            {
                constexpr auto INVALID_INDEX = (size_t) - 1;
                const CM::Point2 ab {anEndPoint.myX - aStartPoint.myX, anEndPoint.myY - aStartPoint.myY, INVALID_INDEX};
                const CM::Point2 ap {aPoint.myX - aStartPoint.myX, aPoint.myY - aStartPoint.myY, INVALID_INDEX};
                const auto dot = CM::Dot(ab, ap);
                const auto result = (dot * dot) / dist2;
                return std::clamp(result, 0.l, dist2);
            }
        }
    }
}

