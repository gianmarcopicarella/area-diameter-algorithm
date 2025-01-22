//
// Created by Gianmarco prevPcarella on 28/12/24.
//

#include "Utils.h"
#include "CustomMath.h"

#include <cassert>
#include <iostream>
#include <functional>

#define MOD(x, n) ((x) % (n))
#define NEXT(x, n) ((x + 1) % (n))
#define NEXT_CW_IDX(idx, seq) (((idx) + 1) % (seq.size()))
#define PREV_CW_IDX(idx, seq) (((idx) - 1) % (seq.size()))

namespace MT
{
    namespace
    {
        int locCountPtsInTri(const std::vector<CM::Point2>& somePoints, const CM::Point2& p1, const CM::Point2& p2, const CM::Point2& p3, int t = 0)
        {
            int count = 0;
            for(const auto& pt : somePoints)
            {
                if(pt.index != p1.index && pt.index != p2.index && pt.index != p3.index)
                {
                    if( CM::AreCollinear(p1, p2, pt) ||
                        CM::AreCollinear(p1, p3, pt) ||
                        CM::AreCollinear(p3, p2, pt))
                    {
                        //count += 1;
                        continue;
                        //continue;
                    }

                    long double s1 = (pt.x - p2.x) * (p1.y - p2.y) - (p1.x - p2.x) * (pt.y - p2.y);
                    long double s2 = (pt.x - p3.x) * (p2.y - p3.y) - (p2.x - p3.x) * (pt.y - p3.y);
                    long double s3 = (pt.x - p1.x) * (p3.y - p1.y) - (p3.x - p1.x) * (pt.y - p1.y);
                    if(CM::IsCloseToZero(s1)) s1=0;
                    if(CM::IsCloseToZero(s2)) s2=0;
                    if(CM::IsCloseToZero(s3)) s3=0;

                    if((s1 < 0.f && s2 < 0.f && s3 < 0.f) ||
                       (s1 > 0.f && s2 > 0.f && s3 > 0.f))
                    {
                        if(t == 1)
                            std::cout << pt.x << ", " << pt.y << ", " << pt.index <<std::endl;
                        ++count;
                    }
                }
            }
            return count;
        }

        int locGetFirstClockWiseLeftIndex(const std::vector<CM::Point2>& somePoints, const CM::Point2& aRefPoint)
        {
            for(int i = 0; i < somePoints.size() - 1; ++i)
            {
                if(somePoints[i].x > aRefPoint.x && somePoints[i+1].x <= aRefPoint.x)
                {
                    return i + 1;
                }
            }
            return 0;
        }
    }

    void CountPointsBelowAllSegments(const std::vector<CM::Point2>& somePoints,
                                     const std::vector<std::vector<CM::Point2>>& someClockWiseSortedPoints,
                                     std::vector<std::vector<int>>& someOutBelowCounts,
                                     std::vector<std::vector<int>>& someOutCollinearCounts)
    {
        // Sort points by increasing x-coordinate
        std::vector<CM::Point2> sortedPoints(somePoints);
        std::sort(sortedPoints.begin(), sortedPoints.end(),CM::SortPointsHorizontally);

        for (int i = 1; i < sortedPoints.size(); ++i)
        {
            const auto& refP = sortedPoints[i];
            const auto& clockwisePoints = someClockWiseSortedPoints[refP.index];
            const auto startIdx = locGetFirstClockWiseLeftIndex(clockwisePoints, refP);

            for(int j = NEXT_CW_IDX(startIdx, clockwisePoints), jp = startIdx;
                j != MOD(startIdx + i, clockwisePoints.size());
                jp = j, j = NEXT_CW_IDX(j, clockwisePoints))
            {
                const auto& prevP = clockwisePoints[jp];
                const auto& currP = clockwisePoints[j];
                const auto areCollinear = CM::AreCollinear(refP, prevP, currP);

                if(areCollinear)
                {
                    someOutCollinearCounts[refP.index][currP.index] = 1 + someOutCollinearCounts[refP.index][prevP.index];
                    someOutCollinearCounts[currP.index][refP.index] = someOutCollinearCounts[refP.index][currP.index];
                }

                if(currP.x < prevP.x && areCollinear)
                {
                    someOutBelowCounts[currP.index][refP.index] = someOutBelowCounts[prevP.index][refP.index] +
                                                                    someOutBelowCounts[currP.index][prevP.index];
                }
                else if(currP.x < prevP.x)
                {
                    const auto b = someOutBelowCounts[prevP.index][refP.index];
                    const auto a = someOutBelowCounts[currP.index][prevP.index];
                    const auto c1 = someOutCollinearCounts[refP.index][prevP.index];
                    const auto c2 = someOutCollinearCounts[prevP.index][currP.index];

                    someOutBelowCounts[currP.index][refP.index] = b + a + c1 + c2 + 1;
                }
                else if(currP.x > prevP.x)
                {
                    someOutBelowCounts[currP.index][refP.index] =
                            (someOutBelowCounts[prevP.index][refP.index] + someOutCollinearCounts[prevP.index][refP.index]) -
                            (someOutBelowCounts[prevP.index][currP.index] + someOutCollinearCounts[prevP.index][currP.index]);
                }
                someOutBelowCounts[refP.index][currP.index] = someOutBelowCounts[currP.index][refP.index];

                // Check the count is correct
                /*
                int countBelow = 0;
                for(auto tp : somePoints)
                {
                    if( tp.index != refP.index &&
                        tp.index != currP.index &&
                        tp.x >= currP.x && tp.x <= refP.x &&
                        CM::Orientation(refP, currP, tp) == CM::ORIENTATION::LEFT_TURN)
                    {
                        std::cout << tp.index << std::endl;
                        ++countBelow;
                    }
                }
                std::cout << "---" << std::endl;

                const auto count1 = someOutBelowCounts[refP.index][currP.index];
                if(count1 != countBelow)
                {
                    const auto a1 = locProcessAngle(refP, prevP);
                    const auto a2 = locProcessAngle(refP, currP);
                    std::cout << "wtf: " << a1 << ", " << a2 << ", " << locClockWiseCompare(refP, prevP, currP) << std::endl;
                }
                assert(countBelow == count1);
                */
            }
        }

/*
        for(int i = 0; i < somePoints.size(); ++i)
        {
            for(int j = 0; j < somePoints.size(); ++j)
            {
                for(int l = 0; l < somePoints.size(); ++l)
                {
                    if(i != j && i != l && j != l)
                    {
                        const auto count1 = locCountPtsInTri(somePoints, somePoints[i], somePoints[j], somePoints[l]);
                        const auto count2 = PointsInTriangle(somePoints[i], somePoints[j], somePoints[l], someOutBelowCounts, someOutCollinearCounts);
                        if(count1 != count2)
                        {
                            locCountPtsInTri(somePoints, somePoints[i], somePoints[j], somePoints[l], 1);
                            PointsInTriangle(somePoints[i], somePoints[j], somePoints[l], someOutBelowCounts, someOutCollinearCounts);
                        }

                        //std::cout << count1 << " and " << count2 << "for tri " << i << ", " << j << ", " << l << std::endl;
                        assert(count1 == count2);
                    }
                }
            }
        }
*/
    }

    int PointsInTriangle(
            const CM::Point2& aFirstPoint,
            const CM::Point2& aSecondPoint,
            const CM::Point2& aThirdPoint,
            const std::vector<std::vector<int>>& somePointCountBelowSegments,
            const std::vector<std::vector<int>>& someCollinearPointCounts)
    {
        if(CM::AreCollinear(aFirstPoint, aSecondPoint, aThirdPoint))
        {
            return 0;
        }

        CM::Point2 points[3] = {aFirstPoint, aSecondPoint, aThirdPoint};
        std::sort(points, points + 3, CM::SortPointsHorizontally);
        const auto li = points[0].index;
        const auto mi = points[1].index;
        const auto ri = points[2].index;
        const auto count = std::abs(
                somePointCountBelowSegments[li][mi] +
                somePointCountBelowSegments[mi][ri] -
                somePointCountBelowSegments[li][ri]);

        const auto orientation = CM::Orientation(points[2], points[0], points[1]);
        if(orientation == CM::ORIENTATION::LEFT_TURN)
        {
            // not counting collinear points
            return count - 1 - someCollinearPointCounts[li][mi] - someCollinearPointCounts[mi][ri];
            // counting collinear points
            // return  count - 1 + someCollinearPointCounts[li][ri];
        }
        // not counting collinear points
        return count - someCollinearPointCounts[li][ri];
        // counting collinear points
        // return count + someCollinearPointCounts[li][mi] + someCollinearPointCounts[mi][ri];
    }

    bool ArePointsClockwise(const CM::Point2& aReferencePoint, const CM::Point2& aFirstPoint, const CM::Point2& aSecondPoint)
    {
        assert(aFirstPoint.index != aSecondPoint.index);
        const auto firstAngle = CM::Angle(aReferencePoint, aFirstPoint);
        const auto secondAngle = CM::Angle(aReferencePoint, aSecondPoint);
        if(CM::IsCloseToZero(firstAngle - secondAngle))
        {
            return CM::SquaredDistance(aReferencePoint, aFirstPoint) <
                   CM::SquaredDistance(aReferencePoint, aSecondPoint);
        }
        return firstAngle > secondAngle;
    }

    void SortPointsClockWiseAroundPoint(const CM::Point2& aReferencePoint, std::vector<CM::Point2>& someOutClockWiseSortedPoints)
    {
        using namespace std::placeholders;
        std::sort(someOutClockWiseSortedPoints.begin(), someOutClockWiseSortedPoints.end(),
                  std::bind(ArePointsClockwise, aReferencePoint, _1, _2));
    }

    long double ComputeDiameter(const std::vector<CM::Point2>& somePoints)
    {
        auto result = 0.l;
        for (const auto& point : somePoints)
        {
            for (const auto& anotherPoint : somePoints)
            {
                if(point.index != anotherPoint.index)
                {
                    const auto distance2 = CM::SquaredDistance(point, anotherPoint);
                    if(distance2 > result)
                    {
                        result = distance2;
                    }
                }
            }
        }
        return std::sqrtl(result);
    }

    void FindAntipodalPairs(
            const std::vector<CM::Point2>& somePoints,
            std::vector<std::pair<size_t, size_t>>& someOutAntipodalIndices)
    {
        assert(somePoints.size() > 2);
        const auto pointsCount = somePoints.size();
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
            someOutAntipodalIndices.emplace_back(somePoints[pi].index, somePoints[qi].index);

            while(  CM::SignedArea(somePoints[pi], somePoints[NEXT(pi, pointsCount)], somePoints[NEXT(qi, pointsCount)]) >
                    CM::SignedArea(somePoints[pi], somePoints[NEXT(pi, pointsCount)], somePoints[qi]))
            {
                qi = NEXT(qi, pointsCount);
                if(pi != q0 || qi != 0)
                {
                    someOutAntipodalIndices.emplace_back(somePoints[pi].index, somePoints[qi].index);
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
                    someOutAntipodalIndices.emplace_back(somePoints[pi].index, somePoints[NEXT(qi, pointsCount)].index);
                }
                else
                {
                    break;
                }
            }
        }
    }
}
