//
// Created by Gianmarco Picarella on 28/12/24.
//

#include "Utils.h"
#include "CustomMath.h"

#include <cassert>
#include <iostream>
#include <functional>

#define MOD(x, n) ((x) % (n))
#define NEXT(x, n) ((x + 1) % (n))
#define NEXT_CW_IDX(idx, seq) (((idx) + 1) % (seq.size()))

//#define CHECK_POINTS_COUNT_CORRECTNESS

namespace MT
{
    namespace
    {
#ifdef CHECK_POINTS_COUNT_CORRECTNESS
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
                        //if(t == 1)
                        //    std::cout << pt.myX << ", " << pt.myY << ", " << pt.myIndex <<std::endl;
                        ++count;
                    }
                }
            }
            return count;
        }
#endif
        int locGetFirstClockWiseLeftIndex(const std::vector<CM::Point2>& somePoints, const CM::Point2& aPoint)
        {
            for(int i = 0; i < somePoints.size() - 1; ++i)
            {
                if(somePoints[i].myX > aPoint.myX && somePoints[i+1].myX <= aPoint.myX)
                {
                    return i + 1;
                }
            }
            return 0;
        }
    }

    bool Diameter::operator==(const Diameter& anotherDiameter) const
    {
        return  myFirstIndex == anotherDiameter.myFirstIndex && mySecondIndex == anotherDiameter.mySecondIndex ||
                mySecondIndex == anotherDiameter.myFirstIndex && myFirstIndex == anotherDiameter.mySecondIndex;
    }

    void CountPointsBelowAllSegments(const std::vector<CM::Point2>& somePoints,
                                     const std::vector<std::vector<CM::Point2>>& someClockWiseSortedPoints,
                                     PointsInTriangleCache& anOutCache)
    {
        anOutCache.myPointsBelowSegmentCount = std::vector<std::vector<int>>(somePoints.size(), std::vector<int>(somePoints.size(), 0));
        anOutCache.myCollinearPointsCount = std::vector<std::vector<int>>(somePoints.size(), std::vector<int>(somePoints.size(), 0));
        anOutCache.myPointsRightBelowCount = std::vector<int>(somePoints.size(), 0);

        // Sort points by increasing x-coordinate
        std::vector<CM::Point2> sortedPoints(somePoints);
        std::sort(sortedPoints.begin(), sortedPoints.end(), CM::SortPointsHorizontally);

        for (int i = 1; i < sortedPoints.size(); ++i)
        {
            const auto& pi = sortedPoints[i];
            const auto& clockwisePoints = someClockWiseSortedPoints[pi.myIndex];
            const auto startIdx = locGetFirstClockWiseLeftIndex(clockwisePoints, pi);

            if(pi.myX == sortedPoints[i - 1].myX)
            {
                anOutCache.myPointsRightBelowCount[pi.myIndex] = 1 + anOutCache.myPointsRightBelowCount[sortedPoints[i - 1].myIndex];
            }

            for(int j = NEXT_CW_IDX(startIdx, clockwisePoints), jp = startIdx;
                j != MOD(startIdx + i, clockwisePoints.size());
                jp = j, j = NEXT_CW_IDX(j, clockwisePoints))
            {
                const auto& pjPrev = clockwisePoints[jp];
                const auto& pj = clockwisePoints[j];

                assert(CM::Orientation(pi, pj, pjPrev) >= CM::ORIENTATION::COLLINEAR);

                if(CM::AreCollinear(pi, pjPrev, pj))
                {
                    anOutCache.myCollinearPointsCount[pi.myIndex][pj.myIndex] = 1 + anOutCache.myCollinearPointsCount[pi.myIndex][pjPrev.myIndex];
                    if(pj.myX == pi.myX) continue;

                    const auto a = anOutCache.myPointsBelowSegmentCount[pi.myIndex][pjPrev.myIndex];
                    const auto b = anOutCache.myPointsBelowSegmentCount[pjPrev.myIndex][pj.myIndex];
                    const auto c = anOutCache.myPointsRightBelowCount[pjPrev.myIndex];
                    anOutCache.myPointsBelowSegmentCount[pj.myIndex][pi.myIndex] = a + b + c;
                }
                else if(pj.myX > pjPrev.myX)
                {
                    if(pj.myX == pi.myX) continue;
                    const auto a = anOutCache.myPointsBelowSegmentCount[pi.myIndex][pjPrev.myIndex];
                    const auto b = anOutCache.myCollinearPointsCount[pi.myIndex][pjPrev.myIndex];
                    const auto c = anOutCache.myPointsBelowSegmentCount[pj.myIndex][pjPrev.myIndex];
                    const auto d = anOutCache.myPointsRightBelowCount[pj.myIndex];

                    anOutCache.myPointsBelowSegmentCount[pj.myIndex][pi.myIndex] = (a + b) - c - d;
                }
                else if(pj.myX < pjPrev.myX)
                {
                    const auto a = anOutCache.myPointsBelowSegmentCount[pi.myIndex][pjPrev.myIndex];
                    const auto b = anOutCache.myCollinearPointsCount[pi.myIndex][pjPrev.myIndex];
                    const auto c = anOutCache.myPointsBelowSegmentCount[pjPrev.myIndex][pj.myIndex];
                    const auto d = anOutCache.myCollinearPointsCount[pjPrev.myIndex][pj.myIndex];
                    const auto e = anOutCache.myPointsRightBelowCount[pjPrev.myIndex];

                    anOutCache.myPointsBelowSegmentCount[pj.myIndex][pi.myIndex] = a + b + c + d + e;
                    if(pjPrev.myX != pi.myX) anOutCache.myPointsBelowSegmentCount[pj.myIndex][pi.myIndex] += 1;
                }
                else
                {
                    anOutCache.myPointsBelowSegmentCount[pj.myIndex][pi.myIndex] = anOutCache.myPointsBelowSegmentCount[pjPrev.myIndex][pi.myIndex];
                }

                anOutCache.myCollinearPointsCount[pj.myIndex][pi.myIndex] = anOutCache.myCollinearPointsCount[pi.myIndex][pj.myIndex];
                anOutCache.myPointsBelowSegmentCount[pi.myIndex][pj.myIndex] = anOutCache.myPointsBelowSegmentCount[pj.myIndex][pi.myIndex];

#ifdef CHECK_POINTS_COUNT_CORRECTNESS
                // Check that collinear counts are correct
                {
                    size_t collinearCount = 0;
                    for(const auto& p : somePoints)
                    {
                        if(p.myIndex != pi.myIndex && p.myIndex != pj.myIndex && p.myX <= pi.myX && CM::AreCollinear(pi, pj, p))
                        {
                            const CM::Point2 first {pj.myX - pi.myX, pj.myY - pi.myY, INVALID_INDEX };
                            const CM::Point2 second {p.myX - pi.myX, p.myY - pi.myY, INVALID_INDEX };
                            const auto dot = CM::Dot(first, second);
                            if(dot > 0 && dot < CM::Distance2(pi, pj))
                            {
                                // std::cout << pi.myIndex << ", " << pj.myIndex << ", " << p.myIndex << ", " << pjPrev.myIndex << ", " << CM::Distance2(pi, p) << ", " << CM::Distance2(pi, pj) << std::endl;
                                ++collinearCount;
                            }
                        }
                    }

                    const auto myCollinearCount = anOutCache.myCollinearPointsCount[pi.myIndex][pj.myIndex];
                    if(myCollinearCount != collinearCount)
                    {
                        std::cout << "wtf: " << collinearCount << ", " << myCollinearCount <<std::endl;
                        std::cout << pi.myIndex << ", " << pj.myIndex << ", " << pjPrev.myIndex << std::endl;
                        throw std::runtime_error("WTF coll!!!");
                    }
                    assert(collinearCount == myCollinearCount);
                }

                // Check that below counts are correct
                {
                    size_t belowCount = 0;
                    for(const auto& p : somePoints)
                    {
                        if( p.myIndex != pi.myIndex &&
                            p.myIndex != pj.myIndex &&
                            p.myX > pj.myX && p.myX < pi.myX &&
                            CM::Orientation(pi, pj, p) == CM::ORIENTATION::LEFT_TURN)
                        {
                            // std::cout << "(counter)" << pi.myIndex << ", " << pj.myIndex << ", " << p.myIndex << std::endl;
                            ++belowCount;
                        }
                    }

                    //std::cout << "----" << std::endl;

                    auto myBelowCount { 0 };
                    std::string branch;
                    if(pi.myX == pj.myX)
                    {
                        if(pi.myY < pj.myY)
                        {
                            myBelowCount =
                                    anOutCache.myPointsRightBelowCount[pj.myIndex] -
                                    (anOutCache.myPointsRightBelowCount[pi.myIndex] + 1);
                            branch = "1";
                        }
                        else
                        {
                            myBelowCount =
                                    anOutCache.myPointsRightBelowCount[pi.myIndex] -
                                    (anOutCache.myPointsRightBelowCount[pj.myIndex] + 1);
                            branch = "2";
                        }
                    }
                    else
                    {
                        myBelowCount = anOutCache.myPointsBelowSegmentCount[pi.myIndex][pj.myIndex];
                        branch = "3";
                    }

                    if(myBelowCount != belowCount)
                    {
                        std::cout << "wtf: " << belowCount << ", " << myBelowCount << " branch=" << branch << std::endl;
                        std::cout << pi.myIndex << ", " << pj.myIndex << ", " << pjPrev.myIndex << ", " << (pi.myX == pjPrev.myX) << std::endl;
                        throw std::runtime_error("WTF below!!!");
                    }
                    assert(belowCount == myBelowCount);
                }
#endif
            }
        }

#ifdef CHECK_POINTS_COUNT_CORRECTNESS
        for(int i = 0; i < somePoints.size(); ++i)
        {
            for(int j = i+1; j < somePoints.size(); ++j)
            {
                for(int l = j+1; l < somePoints.size(); ++l)
                {
                    const auto count1 = locCountPtsInTri(somePoints, somePoints[i], somePoints[j], somePoints[l]);
                    const auto count2 = PointsInTriangle(somePoints[i], somePoints[j], somePoints[l], anOutCache);
                    if(count1 != count2)
                    {
                        std::cout << count1 << ", " << count2 << ", " << i << ", " << j << ", " << l << std::endl;
                        locCountPtsInTri(somePoints, somePoints[i], somePoints[j], somePoints[l], 1);
                        PointsInTriangle(somePoints[i], somePoints[j], somePoints[l], anOutCache);
                        throw std::runtime_error("!!!");
                    }
                    assert(count1 == count2);
                }
            }
        }
#endif
    }

    int PointsInTriangle(
            const CM::Point2& aFirstPoint,
            const CM::Point2& aSecondPoint,
            const CM::Point2& aThirdPoint,
            const PointsInTriangleCache& aCache)
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

        if(points[0].myX == points[1].myX)
        {
            return  aCache.myPointsBelowSegmentCount[mi][ri] -
                    aCache.myPointsBelowSegmentCount[li][ri] -
                    aCache.myCollinearPointsCount[li][ri];
        }
        else if(points[1].myX == points[2].myX)
        {
            return  aCache.myPointsBelowSegmentCount[li][ri] -
                    aCache.myPointsBelowSegmentCount[li][mi] -
                    aCache.myCollinearPointsCount[li][mi];
        }
        else if(CM::Orientation(points[2], points[0], points[1]) == CM::ORIENTATION::LEFT_TURN)
        {
            return  aCache.myPointsBelowSegmentCount[li][ri] - aCache.myPointsBelowSegmentCount[li][mi] -
                    aCache.myPointsBelowSegmentCount[mi][ri] - aCache.myCollinearPointsCount[li][mi] -
                    aCache.myCollinearPointsCount[mi][ri] - aCache.myPointsRightBelowCount[mi] - 1;
        }
        else
        {
            return  aCache.myPointsBelowSegmentCount[li][mi] + aCache.myPointsBelowSegmentCount[mi][ri] +
                    aCache.myPointsRightBelowCount[mi] - aCache.myPointsBelowSegmentCount[li][ri] -
                    aCache.myCollinearPointsCount[li][ri];
        }
    }

    bool ArePointsClockwise(const CM::Point2& aReferencePoint, const CM::Point2& aFirstPoint, const CM::Point2& aSecondPoint)
    {
        assert(aFirstPoint.myIndex != aSecondPoint.myIndex);
        const auto firstAngle = CM::Angle(aReferencePoint, aFirstPoint);
        const auto secondAngle = CM::Angle(aReferencePoint, aSecondPoint);

        if( std::fabsl(firstAngle - secondAngle) < 0.000001l &&
            CM::AreCollinear(aReferencePoint, aFirstPoint, aSecondPoint))
        {
            return  CM::Distance2(aReferencePoint, aFirstPoint) <
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

    std::optional<Diameter> ComputeDiameter(const std::vector<CM::Point2>& somePoints)
    {
        if(somePoints.size() < 2)
        {
            return std::nullopt;
        }

        std::pair<size_t, size_t> bestPair {0, 0};
        long double bestDistance2 { 0 };

        ForAllAntipodalPairs(somePoints, [&](const size_t aFirstIndex, const size_t aSecondIndex){
            const auto distance2 = CM::Distance2(somePoints[aFirstIndex], somePoints[aSecondIndex]);
            if(distance2 > bestDistance2)
            {
                bestDistance2 = distance2;
                bestPair = { aFirstIndex, aSecondIndex };
            }
        });

        return Diameter { somePoints[bestPair.first].myIndex, somePoints[bestPair.second].myIndex };
    }

    void GetHullPoints(const std::vector<size_t>& someIndices,
                       const std::vector<CM::Point2>& somePoints,
                       std::vector<CM::Point2>& someOutHullPoints)
    {
        someOutHullPoints.resize(someIndices.size());
        std::transform(someIndices.begin(), someIndices.end(), someOutHullPoints.begin(),
                       [&](auto index){ return somePoints[index]; });
    }
}

