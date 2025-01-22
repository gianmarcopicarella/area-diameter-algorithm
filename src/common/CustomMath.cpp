//
// Created by Gianmarco Picarella on 27/12/24.
//

#include "CustomMath.h"

#include <cmath>
#include <functional>

namespace MT
{
    namespace CM
    {
        long double SquaredDistance(const Point2& aFirstPoint, const Point2& aSecondPoint)
        {
            const long double dx = aSecondPoint.x - aFirstPoint.x;
            const long double dy = aSecondPoint.y - aFirstPoint.y;
            return dx * dx + dy * dy;
        }

        bool IsCloseToZero(long double aValue, long double anEpsilon)
        {
            return std::fabsl(aValue) < anEpsilon;
        }

        long double Dot2(const Point2& aFirstPoint, const Point2& aSecondPoint)
        {
            return aFirstPoint.x * aSecondPoint.x + aFirstPoint.y * aSecondPoint.y;
        }

        long double SignedArea(const Point2& aFirstPoint, const Point2& aSecondPoint, const Point2& aThirdPoint)
        {
            return 0.5l * ( (aSecondPoint.x - aFirstPoint.x) * (aThirdPoint.y - aFirstPoint.y) -
                            (aThirdPoint.x - aFirstPoint.x) * (aSecondPoint.y - aFirstPoint.y));
        }

        bool AreCollinear(const Point2& aFirstPoint, const Point2& aSecondPoint, const Point2& aThirdPoint)
        {
            return IsCloseToZero(SignedArea(aFirstPoint, aSecondPoint, aThirdPoint));
        }

        ORIENTATION Orientation(const Point2& aFirstPoint, const Point2& aSecondPoint, const Point2& aThirdPoint)
        {
            const auto signedArea = SignedArea(aFirstPoint, aSecondPoint, aThirdPoint);
            if(IsCloseToZero(signedArea))
            {
                return ORIENTATION::COLLINEAR;
            }
            else if(signedArea > 0)
            {
                return ORIENTATION::LEFT_TURN;
            }
            else
            {
                return ORIENTATION::RIGHT_TURN;
            }
        }

        long double Angle(const CM::Point2& aReferencePoint, const CM::Point2& aPoint)
        {
            const auto angle =
                    std::atan2l(aPoint.y - aReferencePoint.y, aPoint.x - aReferencePoint.x);
            if(angle < 0.f)
            {
                return angle + 2.f * M_PI;
            }
            return angle;
        }
    }
}

