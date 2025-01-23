//
// Created by Gianmarco Picarella on 06/01/25.
//

#include <iostream>

#include "Antipodal.h"
#include "Utils.h"

// https://tessil.github.io/2016/08/29/benchmark-hopscotch-map.html
// #define TSL_NO_EXCEPTIONS
#include <tsl/hopscotch_map.h>
#include <tsl/bhopscotch_map.h>

#define KEY(a, b, c, d, n) ((a) * (n) * (n) * (n) + (b) * (n) * (n) + (c) * (n) + (d))
// #define MERGE(a, b, c, d) (((uint64_t)a) + (((uint64_t)b) << 16) + (((uint64_t)c) << 32) + (((uint64_t)d) << 48));
//#define POINT_FILTER_TRI_AREA
//#define

namespace MT
{
    namespace
    {
        void locPartitionLeftAndRightPoints(
                const std::vector<CM::Point2>& somePoints,
                const CM::Point2& aStartPoint,
                const CM::Point2& anEndPoint,
                std::vector<CM::Point2>& someOutLeftPoints,
                std::vector<CM::Point2>& someOutRightPoints)
        {
            const auto maximumDistance2 = CM::SquaredDistance(aStartPoint, anEndPoint);
            constexpr auto INVALID_INDEX = (size_t) - 1;
            const CM::Point2 startToEnd { anEndPoint.x - aStartPoint.x, anEndPoint.y - aStartPoint.y, INVALID_INDEX };
            const CM::Point2 endToStart { aStartPoint.x - anEndPoint.x, aStartPoint.y - anEndPoint.y, INVALID_INDEX };
            for (const auto& point : somePoints)
            {
                if( point.index != aStartPoint.index &&
                    point.index != anEndPoint.index &&
                    CM::SquaredDistance(aStartPoint, point) <= maximumDistance2 &&
                    CM::SquaredDistance(anEndPoint, point) <= maximumDistance2)
                {
                    const CM::Point2 startToPoint { point.x - aStartPoint.x, point.y - aStartPoint.y, INVALID_INDEX };
                    const CM::Point2 endToPoint { point.x - anEndPoint.x, point.y - anEndPoint.y, INVALID_INDEX };
                    if( CM::Dot2(startToEnd, startToPoint) >= 0 &&
                        CM::Dot2(endToStart, endToPoint) >= 0)
                    {
                        if(CM::Orientation(aStartPoint, anEndPoint, point) <= CM::ORIENTATION::COLLINEAR)
                        {
                            someOutRightPoints.emplace_back(point);
                        }
                        else
                        {
                            someOutLeftPoints.emplace_back(point);
                        }
                    }
                }
            }
        }

        void locKeepUniquePolyPoints(
            const CM::Point2& s, const CM::Point2& t,
            const CM::Point2& l, const CM::Point2& pl,
            const CM::Point2& r, const CM::Point2& pr,
            std::array<CM::Point2, 6>& someOutPoints,
            size_t& anOutLength)
        {
            anOutLength = 0;
            someOutPoints[anOutLength++] = t;
            if(pl.index != someOutPoints[anOutLength - 1].index)
            {
                someOutPoints[anOutLength++] = pl;
            }
            if(l.index != someOutPoints[anOutLength - 1].index)
            {
                someOutPoints[anOutLength++] = l;
            }
            if(s.index != someOutPoints[anOutLength - 1].index)
            {
                someOutPoints[anOutLength++] = s;
            }
            if(pr.index != someOutPoints[anOutLength - 1].index)
            {
                someOutPoints[anOutLength++] = pr;
            }
            if( r.index != someOutPoints[anOutLength - 1].index &&
                r.index != someOutPoints[0].index)
            {
                someOutPoints[anOutLength++] = r;
            }
        }

        bool locAreAntipodal(const std::vector<std::pair<size_t, size_t>>& somePairs,
                             const CM::Point2& aLeftPoint, const CM::Point2& aRightPoint)
        {
            return std::find_if(somePairs.begin(), somePairs.end(), [&](const auto& pair){
                return  (pair.first == aLeftPoint.index && pair.second == aRightPoint.index) ||
                        (pair.first == aRightPoint.index && pair.second == aLeftPoint.index);
            }) != somePairs.end();
        }

        bool locIsValidEdge(const CM::Point2& a, const CM::Point2& b, const CM::Point2& c, const CM::Point2& s, const CM::Point2& t)
        {
            return  CM::Orientation(a, b, s) >= CM::ORIENTATION::COLLINEAR &&
                    CM::Orientation(a, b, t) >= CM::ORIENTATION::COLLINEAR &&
                    CM::Orientation(a, b, c) >= CM::ORIENTATION::COLLINEAR;
        }

        void locCopyListAndPoint(const std::vector<CM::Point2>& somePointsToCopy,
                                 const CM::Point2& anotherPoint,
                                 std::vector<CM::Point2>& someOutPoints)
        {
            someOutPoints.resize(somePointsToCopy.size() + 1);
            std::copy(somePointsToCopy.begin(), somePointsToCopy.end(), someOutPoints.begin());
            someOutPoints.back() = anotherPoint;
        }

        void locInsertOrUpdateMinArea(const size_t aCacheKey,
                                      const long double aCurrentArea,
                                      const long double aMaxArea,
                                      tsl::hopscotch_map<size_t, long double>& anOutcache)
        {
            if(aCurrentArea <= aMaxArea)
            {
                auto nextAreaIter = anOutcache.find(aCacheKey);
                if(nextAreaIter == anOutcache.end())
                {
                    anOutcache.insert({aCacheKey, aCurrentArea});
                }
                else
                {
                    nextAreaIter.value() = std::fminl(nextAreaIter.value(), aCurrentArea);
                }
            }
        }

        AntipodalResult locProcessSegment(const CM::Point2& aStartPoint,
                                        const CM::Point2& anEndPoint,
                                        const std::vector<CM::Point2>& somePoints,
                                        const std::vector<CM::Point2>& someLeftPoints,
                                        const std::vector<CM::Point2>& someRightPoints,
                                        const std::vector<std::vector<int>>& someBelowPointsCounts,
                                        const std::vector<std::vector<int>>& someCollinearPointsCounts,
                                        const size_t aTotalPointsCount,
                                        const size_t aMaxPointsCount,
                                        const long double aMaxArea)
        {
            const auto maxAllowedPoints = std::min(aMaxPointsCount, someRightPoints.size() + someLeftPoints.size());
            if(maxAllowedPoints == 0) return AntipodalResult{ std::numeric_limits<long double>::infinity(), 0 };
            // tsl::hopscotch_map<size_t, long double> cache;
            // std::vector<std::vector<std::tuple<CM::Point2, CM::Point2, CM::Point2, CM::Point2>>> iterCache {maxAllowedPoints + 1, std::vector<std::tuple<CM::Point2, CM::Point2, CM::Point2, CM::Point2>>{}};

            struct Entry { size_t l, pl, r, pr, q; long double a; };
            std::vector<std::unordered_map<size_t, Entry>> cache {maxAllowedPoints + 1,
                                                                  std::unordered_map<size_t, Entry>{}};

            std::vector<CM::Point2> leftPointsWithStart, leftPointsWithEnd,
                                    rightPointsWithStart, rightPointsWithEnd;

            locCopyListAndPoint(someLeftPoints, aStartPoint, leftPointsWithStart);
            locCopyListAndPoint(someRightPoints, aStartPoint, rightPointsWithStart);
            locCopyListAndPoint(someLeftPoints, anEndPoint, leftPointsWithEnd);
            locCopyListAndPoint(someRightPoints, anEndPoint, rightPointsWithEnd);

            // Initialize entries at k = 0
            for (const auto& pl: leftPointsWithEnd)
            {
                for (const auto& pr: rightPointsWithStart)
                {
                    const auto key = KEY(aStartPoint.index, pl.index, anEndPoint.index, pr.index, aTotalPointsCount);
                    Entry entry {aStartPoint.index, pl.index, anEndPoint.index, pr.index, 0, 0};
                    cache[0].insert({key, entry});
                }
            }

            const auto segmentLength2 = CM::SquaredDistance(aStartPoint, anEndPoint);
            AntipodalResult result {std::numeric_limits<long double>::infinity(), 0};


            std::array<CM::Point2, 6> poly;

            for (size_t k = 0; k < maxAllowedPoints + 1; ++k)
            {
                for (const auto& pair : cache[k])
                {
                    const auto& entry = pair.second;
                    const auto& l = somePoints[entry.l];
                    const auto& r = somePoints[entry.r];

                    if(CM::SquaredDistance(l, r) > segmentLength2)
                    {
                        continue;
                    }

                    const auto& pl = somePoints[entry.pl];
                    const auto& pr = somePoints[entry.pr];
                    const auto& currentArea = entry.a;

                    if(pl.index == anEndPoint.index && pr.index == aStartPoint.index)
                    {
                        if( k > 0 &&
                            (result.myPointsCount < k ||
                             (result.myPointsCount == k && currentArea < result.myHullArea)))
                        {
                            result.myHullArea = currentArea;
                            result.myPointsCount = k;
                        }
                        continue;
                    }

                    size_t aPolyLength = 0;
                    locKeepUniquePolyPoints(aStartPoint, anEndPoint, l, pl, r, pr, poly, aPolyLength);

                    /*
                    Should not be needed because the previous check already filters out all the cases for which currPoly
                    has less than 3 unique vertices
                    if(currPoly.size() < 3)
                    {
                        std::cout << "ao" << std::endl;
                        assert(false);
                        continue;
                    }
                    */

                    //std::vector<std::pair<size_t, size_t>> currAntipodalPairs;
                    //FindAntipodalPairs(currPoly, currAntipodalPairs);

                    /*
                    Should not be needed
                    if(!locAreAntipodal(currAntipodalPairs, l, r))
                    {
                        throw std::runtime_error("DIOP!");
                        continue;
                    }
                    */

                    bool isPrevLeftAPWithRight = false, isPrevRightAPWithLeft = false;
                    AreAntipodalPairs(poly, aPolyLength, pl.index, r.index,
                                      pr.index, l.index, isPrevLeftAPWithRight, isPrevRightAPWithLeft);

                    if(isPrevLeftAPWithRight || (pr.index == aStartPoint.index && isPrevRightAPWithLeft))
                    {
                        const auto triangleArea = std::fabsl(CM::SignedArea(anEndPoint, l, pl));
                        const auto triangleCount = 1 + PointsInTriangle(anEndPoint, l, pl, someBelowPointsCounts, someCollinearPointsCounts);

                        for (const auto& ql : leftPointsWithEnd)
                        {
                            if (ql.index != l.index && ql.index != pl.index && locIsValidEdge(ql, pl, l, aStartPoint, anEndPoint))
                            {
                                const auto nextKey = KEY(pl.index, ql.index, r.index, pr.index, aTotalPointsCount);
                                if(currentArea + triangleArea <= aMaxArea)
                                {
                                    auto& nextCache = cache[k + triangleCount];
                                    auto nextAreaPair = nextCache.find(nextKey);
                                    if(nextAreaPair == nextCache.end())
                                    {
                                        Entry newEntry {pl.index, ql.index, r.index, pr.index, l.index, currentArea + triangleArea};
                                        nextCache.insert({nextKey, newEntry});
                                    }
                                    else if(currentArea + triangleArea < nextAreaPair->second.a)
                                    {
                                        nextAreaPair->second.a = currentArea + triangleArea;
                                        nextAreaPair->second.q = l.index;
                                    }
                                }
                            }
                        }
                    }

                    if(isPrevRightAPWithLeft || (pl.index == anEndPoint.index && isPrevLeftAPWithRight))
                    {
                        const auto triangleArea = std::fabsl(CM::SignedArea(aStartPoint, r, pr));
                        const auto triangleCount = 1 + PointsInTriangle(aStartPoint, r, pr, someBelowPointsCounts, someCollinearPointsCounts);

                        for (const auto& qr : rightPointsWithStart)
                        {
                            if (qr.index != r.index && qr.index != pr.index && locIsValidEdge(qr, pr, r, aStartPoint, anEndPoint))
                            {
                                if(currentArea + triangleArea <= aMaxArea)
                                {
                                    const auto nextKey = KEY(l.index, pl.index, pr.index, qr.index, aTotalPointsCount);
                                    if(currentArea + triangleArea <= aMaxArea)
                                    {
                                        auto& nextCache = cache[k + triangleCount];
                                        auto nextAreaPair = nextCache.find(nextKey);
                                        if(nextAreaPair == nextCache.end())
                                        {
                                            Entry newEntry {l.index, pl.index, pr.index, qr.index, r.index, currentArea + triangleArea};
                                            nextCache.insert({nextKey, newEntry});
                                        }
                                        else if(currentArea + triangleArea < nextAreaPair->second.a)
                                        {
                                            nextAreaPair->second.a = currentArea + triangleArea;
                                            nextAreaPair->second.q = r.index;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                /*
                for (const auto& l: leftPointsWithStart)
                {
                    for (const auto& pl: leftPointsWithEnd)
                    {
                        for (const auto& r: rightPointsWithEnd)
                        {
                            for (const auto& pr: rightPointsWithStart)
                            {
                                const auto currKey = TO_KEY(l.index, pl.index, r.index, pr.index, k, aTotalPointsCount);
                                if( CM::SquaredDistance(l, r) > segmentLength2 ||
                                    !cache.contains(currKey))
                                {
                                    continue;
                                }

                                const auto currArea = cache.at(currKey);
                                if(pl.index == anEndPoint.index && pr.index == aStartPoint.index)
                                {
                                    if( k > 0 &&
                                        (result.myPointsCount < k ||
                                        (result.myPointsCount == k && currArea < result.myHullArea)))
                                    {
                                        result.myHullArea = currArea;
                                        result.myPointsCount = k;
                                    }
                                    continue;
                                }

                                std::vector<CM::Point2> currPoly;
                                locKeepUniquePolyPoints(aStartPoint, anEndPoint, l, pl, r, pr, currPoly);
                                if(currPoly.size() < 3)
                                {
                                    continue;
                                }

                                std::vector<std::pair<size_t, size_t>> currAntipodalPairs;
                                FindAntipodalPairs(currPoly, currAntipodalPairs);

                                if(!locAreAntipodal(currAntipodalPairs, l, r))
                                {
                                    continue;
                                }

                                const auto isPrevLeftAPWithRight = locAreAntipodal(currAntipodalPairs, pl, r);
                                const auto isPrevRightAPWithLeft = locAreAntipodal(currAntipodalPairs, pr, l);

                                if(isPrevLeftAPWithRight || (pr.index == aStartPoint.index && isPrevRightAPWithLeft))
                                {
                                    const auto triangleArea = std::fabsl(CM::SignedArea(anEndPoint, l, pl));
                                    const auto triangleCount = 1 + PointsInTriangle(anEndPoint, l, pl, someBelowPointsCounts, someCollinearPointsCounts);
                                    for (const auto& ql : leftPointsWithEnd)
                                    {
                                        if (ql.index != l.index && ql.index != pl.index && locIsValidEdge(ql, pl, l, aStartPoint, anEndPoint))
                                        {
                                            const auto nextKey = TO_KEY(pl.index, ql.index, r.index, pr.index, k + triangleCount, aTotalPointsCount);
                                            locInsertOrUpdateMinArea(nextKey, currArea + triangleArea, aMaxArea, cache);
                                        }
                                    }
                                }

                                if(isPrevRightAPWithLeft || (pl.index == anEndPoint.index && isPrevLeftAPWithRight))
                                {
                                    const auto triangleArea = std::fabsl(CM::SignedArea(aStartPoint, r, pr));
                                    const auto triangleCount = 1 + PointsInTriangle(aStartPoint, r, pr, someBelowPointsCounts, someCollinearPointsCounts);
                                    for (const auto& qr : rightPointsWithStart)
                                    {
                                        if (qr.index != r.index && qr.index != pr.index && locIsValidEdge(qr, pr, r, aStartPoint, anEndPoint))
                                        {
                                            if(currArea + triangleArea <= aMaxArea)
                                            {
                                                const auto nextKey = TO_KEY(l.index, pl.index, pr.index, qr.index, k + triangleCount, aTotalPointsCount);
                                                locInsertOrUpdateMinArea(nextKey, currArea + triangleArea, aMaxArea, cache);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }

                }*/
            }
            return result;
        }
    }

    AntipodalResult AntipodalAlgorithm(
            const std::vector<CM::Point2>& somePoints,
            size_t aMaxPointsCount,
            long double aMaxArea,
            long double aMaxDiameter,
            bool aShouldReconstructHull)
    {
        // Create a 2-dimensional array containing the number of points right below any segment
        const size_t pointsCount = somePoints.size();
        std::vector<std::vector<int>> pointsBelowCounts(pointsCount, std::vector<int>(pointsCount, 0));
        std::vector<std::vector<int>> collinearPointsCounts(pointsCount, std::vector<int>(pointsCount, 0));

        {
            // Create a 2-dimensional array containing for each point the sequence of clockwise sorted points around it
            std::vector<std::vector<CM::Point2>> clockwiseSortedPoints(pointsCount, std::vector<CM::Point2>(pointsCount - 1));
            {
                for(size_t i = 0; i < pointsCount; ++i)
                {
                    std::copy(somePoints.begin(), somePoints.begin() + i, clockwiseSortedPoints[i].begin());
                    std::copy(somePoints.begin() + i + 1, somePoints.end(), clockwiseSortedPoints[i].begin() + i);
                    SortPointsClockWiseAroundPoint(somePoints[i], clockwiseSortedPoints[i]);
                }
            }
            CountPointsBelowAllSegments(somePoints, clockwiseSortedPoints, pointsBelowCounts, collinearPointsCounts);
        }

        const auto maxDiameter2 = aMaxDiameter * aMaxDiameter;
        AntipodalResult result { std::numeric_limits<long double>::infinity(), 0 };

        for (size_t i = 0; i < somePoints.size(); ++i)
        {
            for (size_t j = i + 1; j < somePoints.size(); ++j)
            {
                AntipodalResult pairResult { std::numeric_limits<long double>::infinity(), 0 };
                if(CM::SquaredDistance(somePoints[i], somePoints[j]) <= maxDiameter2)
                {
                    std::vector<CM::Point2> leftPoints, rightPoints;
                    locPartitionLeftAndRightPoints(somePoints,somePoints[i],somePoints[j],leftPoints,rightPoints);

                    /*
                    if(leftPoints.size() + rightPoints.size() < result.myPointsCount)
                    {
                        continue;
                    }
                    */

                    pairResult = locProcessSegment(somePoints[i], somePoints[j], somePoints, leftPoints, rightPoints, pointsBelowCounts, collinearPointsCounts, somePoints.size(), aMaxPointsCount, aMaxArea);
                    if( pairResult.myPointsCount > result.myPointsCount ||
                        (pairResult.myPointsCount == result.myPointsCount && pairResult.myHullArea < result.myHullArea))
                    {
                        result.myPointsCount = pairResult.myPointsCount;
                        result.myHullArea = pairResult.myHullArea;
                    }
                }
#ifdef DEBUG_ANTIPODAL
                if(pairResult.myPointsCount > 0) pairResult.myPointsCount += 2;
                result.results.emplace_back(pairResult.myHullArea, pairResult.myPointsCount);
#endif

#ifdef VERBOSE_ANTIPODAL
                std::cout   << "Done pair (" + std::to_string(i) + ", " + std::to_string(j) + "). Area: " <<
                                pairResult.myHullArea << ", Count: " << pairResult.myPointsCount << std::endl;
#endif
            }
        }

        result.myHasFoundSolution = result.myPointsCount > 0;
        if(result.myHasFoundSolution)
        {
            result.myPointsCount += 2;
            if(aShouldReconstructHull)
            {

            }
        }

        return result;
    }
}