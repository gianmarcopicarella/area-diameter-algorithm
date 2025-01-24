//
// Created by Gianmarco Picarella on 27/12/24.
//

#include "CustomMath.h"

#include <cmath>
#include <cassert>
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
                    std::atan2l(aPoint.y - aReferencePoint.y, aPoint.x - aReferencePoint.x);
            if(angle < 0.f)
            {
                return angle + 2.f * M_PI;
            }
            return angle;
        }

        long double ProjectedLen2(const CM::Point2& aStartPoint, const CM::Point2& anEndPoint, const CM::Point2& aPoint)
        {
            if(aPoint.index == aStartPoint.index)
            {
                return 0;
            }
            const auto dist2 = CM::SquaredDistance(aStartPoint, anEndPoint);
            if(aPoint.index == anEndPoint.index)
            {
                return dist2;
            }
            else
            {
                constexpr auto INVALID_INDEX = (size_t) - 1;
                const CM::Point2 ab {anEndPoint.x - aStartPoint.x, anEndPoint.y - aStartPoint.y, INVALID_INDEX};
                const CM::Point2 ap {aPoint.x - aStartPoint.x, aPoint.y - aStartPoint.y, INVALID_INDEX};
                const auto dot = CM::Dot2(ab, ap);
                const auto result = (dot * dot) / dist2;
                return std::clamp(result, 0.l, dist2);
            }
        }
    }
}

