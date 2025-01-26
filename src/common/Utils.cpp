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
                if(pt.myIndex != p1.myIndex && pt.myIndex != p2.myIndex && pt.myIndex != p3.myIndex)
                {
                    if( CM::AreCollinear(p1, p2, pt) ||
                        CM::AreCollinear(p1, p3, pt) ||
                        CM::AreCollinear(p3, p2, pt))
                    {
                        //count += 1;
                        continue;
                        //continue;
                    }

                    long double s1 = (pt.myX - p2.myX) * (p1.myY - p2.myY) - (p1.myX - p2.myX) * (pt.myY - p2.myY);
                    long double s2 = (pt.myX - p3.myX) * (p2.myY - p3.myY) - (p2.myX - p3.myX) * (pt.myY - p3.myY);
                    long double s3 = (pt.myX - p1.myX) * (p3.myY - p1.myY) - (p3.myX - p1.myX) * (pt.myY - p1.myY);
                    if(CM::IsCloseToZero(s1)) s1=0;
                    if(CM::IsCloseToZero(s2)) s2=0;
                    if(CM::IsCloseToZero(s3)) s3=0;

                    if((s1 < 0.f && s2 < 0.f && s3 < 0.f) ||
                       (s1 > 0.f && s2 > 0.f && s3 > 0.f))
                    {
                        if(t == 1)
                            std::cout << pt.myX << ", " << pt.myY << ", " << pt.myIndex <<std::endl;
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
                if(somePoints[i].myX > aRefPoint.myX && somePoints[i+1].myX <= aRefPoint.myX)
                {
                    return i + 1;
                }
            }
            return 0;
        }
    }


    long double Diameter::Value2() const
    {
        return CM::Distance2(myFirstPoint, mySecondPoint);
    }

    bool Diameter::operator==(const Diameter& anotherDiameter) const
    {
        return  (myFirstPoint == anotherDiameter.myFirstPoint || myFirstPoint == anotherDiameter.mySecondPoint) &&
                (mySecondPoint == anotherDiameter.myFirstPoint || mySecondPoint == anotherDiameter.mySecondPoint);
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
            const auto& clockwisePoints = someClockWiseSortedPoints[refP.myIndex];
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
                    someOutCollinearCounts[refP.myIndex][currP.myIndex] = 1 + someOutCollinearCounts[refP.myIndex][prevP.myIndex];
                    someOutCollinearCounts[currP.myIndex][refP.myIndex] = someOutCollinearCounts[refP.myIndex][currP.myIndex];
                }

                if(currP.myX < prevP.myX && areCollinear)
                {
                    someOutBelowCounts[currP.myIndex][refP.myIndex] = someOutBelowCounts[prevP.myIndex][refP.myIndex] +
                                                                      someOutBelowCounts[currP.myIndex][prevP.myIndex];
                }
                else if(currP.myX < prevP.myX)
                {
                    const auto b = someOutBelowCounts[prevP.myIndex][refP.myIndex];
                    const auto a = someOutBelowCounts[currP.myIndex][prevP.myIndex];
                    const auto c1 = someOutCollinearCounts[refP.myIndex][prevP.myIndex];
                    const auto c2 = someOutCollinearCounts[prevP.myIndex][currP.myIndex];

                    someOutBelowCounts[currP.myIndex][refP.myIndex] = b + a + c1 + c2 + 1;
                }
                else if(currP.myX > prevP.myX)
                {
                    someOutBelowCounts[currP.myIndex][refP.myIndex] =
                            (someOutBelowCounts[prevP.myIndex][refP.myIndex] + someOutCollinearCounts[prevP.myIndex][refP.myIndex]) -
                            (someOutBelowCounts[prevP.myIndex][currP.myIndex] + someOutCollinearCounts[prevP.myIndex][currP.myIndex]);
                }
                someOutBelowCounts[refP.myIndex][currP.myIndex] = someOutBelowCounts[currP.myIndex][refP.myIndex];

                // Check the count is correct
                /*
                int countBelow = 0;
                for(auto tp : somePoints)
                {
                    if( tp.myIndex != refP.myIndex &&
                        tp.myIndex != currP.myIndex &&
                        tp.myX >= currP.myX && tp.myX <= refP.myX &&
                        CM::Orientation(refP, currP, tp) == CM::ORIENTATION::LEFT_TURN)
                    {
                        std::cout << tp.myIndex << std::endl;
                        ++countBelow;
                    }
                }
                std::cout << "---" << std::endl;

                const auto count1 = someOutBelowCounts[refP.myIndex][currP.myIndex];
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
        const auto li = points[0].myIndex;
        const auto mi = points[1].myIndex;
        const auto ri = points[2].myIndex;
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
        assert(aFirstPoint.myIndex != aSecondPoint.myIndex);
        const auto firstAngle = CM::Angle(aReferencePoint, aFirstPoint);
        const auto secondAngle = CM::Angle(aReferencePoint, aSecondPoint);
        if(CM::IsCloseToZero(firstAngle - secondAngle))
        {
            return CM::Distance2(aReferencePoint, aFirstPoint) <
                   CM::Distance2(aReferencePoint, aSecondPoint);
        }
        return firstAngle > secondAngle;
    }

    void SortPointsClockWiseAroundPoint(const CM::Point2& aReferencePoint, std::vector<CM::Point2>& someOutClockWiseSortedPoints)
    {
        using namespace std::placeholders;
        std::sort(someOutClockWiseSortedPoints.begin(), someOutClockWiseSortedPoints.end(),
                  std::bind(ArePointsClockwise, aReferencePoint, _1, _2));
    }

    Diameter ComputeDiameter(const std::vector<CM::Point2>& somePoints)
    {
        Diameter result;
        ForAllAntipodalPairs(somePoints, [&](const auto& p1, const auto& p2){
            if(CM::Distance2(p1, p2) > result.Value2())
            {
                result = {p1, p2};
            }
        });
        return result;
    }
}

