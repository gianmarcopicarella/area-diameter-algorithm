//
// Created by Gianmarco Picarella on 06/01/25.
//

#include <iostream>
#include <queue>

#include "Antipodal.h"
#include "Utils.h"

// https://tessil.github.io/2016/08/29/benchmark-hopscotch-map.html
// #define TSL_NO_EXCEPTIONS
#include <tsl/hopscotch_map.h>
#include <tsl/bhopscotch_map.h>


// #define TO_KEY(a, b, c, d, e, n) ((a) * (n) * (n) * (n) * (n) + (b) * (n) * (n) * (n) + (c) * (n) * (n) + (d) * (n) + (e))
#define KEY(a, b, c, d, n) ((a) * (n) * (n) * (n) + (b) * (n) * (n) + (c) * (n) + (d))
#define KEY2(a, b, c, d, e, n) ((a) * (n) * (n) * (n) * (n) + (b) * (n) * (n) * (n) + (c) * (n) * (n) + (d) * (n) + (e))
//#define POINT_FILTER_TRI_AREA
//#define

namespace MT
{
    namespace
    {
        constexpr auto INVALID_INDEX = (size_t) - 1;

        void locPartitionLeftAndRightPoints(
                const std::vector<CM::Point2>& somePoints,
                const CM::Point2& aStartPoint,
                const CM::Point2& anEndPoint,
                std::vector<CM::Point2>& someOutLeftPoints,
                std::vector<CM::Point2>& someOutRightPoints)
        {
            const auto maximumDistance2 = CM::SquaredDistance(aStartPoint, anEndPoint);
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

            someOutPoints[anOutLength++] = l;

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

        struct Entry { size_t l, pl, r, pr, q; long double a; bool isLeft {false}; };

        template <bool isLeft = false>
        void locInsertOrUpdateMinArea(const size_t a, const size_t b, const size_t c, const size_t d, const size_t e,
                                      const long double aCurrentArea,
                                      const long double aMaxArea,
                                      const size_t aTotalPointsCount,
                                      std::unordered_map<size_t, Entry>& anOutcache)
        {
            if(aCurrentArea <= aMaxArea)
            {
                const auto key = KEY(a, b, c, d, aTotalPointsCount);
                auto nextAreaIter = anOutcache.find(key);
                if(nextAreaIter == anOutcache.end())
                {
                    anOutcache.emplace(key, Entry{a, b, c, d, e, aCurrentArea, isLeft});
                }
                else if(aCurrentArea < nextAreaIter->second.a)
                {
                    nextAreaIter->second.a = aCurrentArea;
                    nextAreaIter->second.q = e;
                    nextAreaIter->second.isLeft = isLeft;
                }
            }
        }

        template <bool Ascending=true>
        void locSortPointsByProjectionLength(
                const CM::Point2& aStartPoint,
                const CM::Point2& anEndPoint,
                std::vector<CM::Point2>& someOutSortedPoints)
        {
            std::sort(someOutSortedPoints.begin(), someOutSortedPoints.end(), [&](auto& a, auto& b){
                if constexpr (Ascending)
                {
                    return
                        CM::ProjectedLen2(aStartPoint, anEndPoint, a) <=
                        CM::ProjectedLen2(aStartPoint, anEndPoint, b);
                }
                else
                {
                    return
                        CM::ProjectedLen2(aStartPoint, anEndPoint, a) >=
                        CM::ProjectedLen2(aStartPoint, anEndPoint, b);
                }
            });
        }

        AntipodalResult locProcessSegment(  const CM::Point2& aStartPoint,
                                            const CM::Point2& anEndPoint,
                                            const std::vector<CM::Point2>& somePoints,
                                            const std::vector<CM::Point2>& someLeftPoints,
                                            const std::vector<CM::Point2>& someRightPoints,
                                            const std::vector<std::vector<int>>& someBelowPointsCounts,
                                            const std::vector<std::vector<int>>& someCollinearPointsCounts,
                                            const size_t aTotalPointsCount,
                                            const size_t aMaxPointsCount,
                                            const long double aMaxArea,
                                            std::vector<std::unordered_map<size_t, Entry>>& cache)
        {
            AntipodalResult result;
            const auto maxAllowedPoints = std::min(aMaxPointsCount, someRightPoints.size() + someLeftPoints.size());

            if(maxAllowedPoints == 0)
            {
                return result;
            }

            std::vector<CM::Point2> leftPointsWithEnd, rightPointsWithStart;
            locCopyListAndPoint(someRightPoints, aStartPoint, rightPointsWithStart);
            locCopyListAndPoint(someLeftPoints, anEndPoint, leftPointsWithEnd);

            // constexpr auto ascending = false;
            std::unordered_map<size_t, size_t> startIndices;
            locSortPointsByProjectionLength<true>(aStartPoint, anEndPoint, leftPointsWithEnd);
            locSortPointsByProjectionLength<false>(aStartPoint, anEndPoint, rightPointsWithStart);
            for(size_t i = 0; i < leftPointsWithEnd.size(); ++i)
            {
                startIndices.insert({leftPointsWithEnd[i].index, i + 1});
            }
            for(size_t i = 0; i < rightPointsWithStart.size(); ++i)
            {
                startIndices.insert({rightPointsWithStart[i].index, i + 1});
            }

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
            std::array<CM::Point2, 6> currentPoly;
            size_t currentPolyLength = 0;

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
                    const auto& currArea = entry.a;

                    if(pl.index == anEndPoint.index && pr.index == aStartPoint.index)
                    {
                        if( k > 0 &&
                            (k > result.myPointsCount ||
                             (result.myPointsCount == k && currArea < result.myHullArea)))
                        {
                            result.myHullArea = currArea;
                            result.myPointsCount = k;
                        }
                        continue;
                    }

                    locKeepUniquePolyPoints(aStartPoint, anEndPoint, l, pl, r, pr, currentPoly, currentPolyLength);

                    bool isPrevLeftAPWithRight = false, isPrevRightAPWithLeft = false;
                    AreAntipodalPairs(currentPoly, currentPolyLength, pl.index,
                                      r.index, l.index, pr.index,
                                      isPrevLeftAPWithRight, isPrevRightAPWithLeft);

                    if(isPrevLeftAPWithRight || (pr.index == aStartPoint.index && isPrevRightAPWithLeft))
                    {
                        const auto triangleArea = std::fabsl(CM::SignedArea(anEndPoint, l, pl));
                        const auto triangleCount = 1 + PointsInTriangle(anEndPoint, l, pl, someBelowPointsCounts, someCollinearPointsCounts);
                        for(size_t i = startIndices[pl.index]; i < leftPointsWithEnd.size(); ++i)
                        {
                            const auto& ql = leftPointsWithEnd[i];
                            assert(ql.index != l.index && ql.index != pl.index);
                            if( CM::Orientation(ql, pl, l) >= CM::ORIENTATION::COLLINEAR &&
                                CM::Orientation(ql, pl, anEndPoint) >= CM::ORIENTATION::COLLINEAR)
                            {
                                locInsertOrUpdateMinArea<true>(pl.index, ql.index, r.index, pr.index, l.index,
                                                               currArea + triangleArea, aMaxArea, aTotalPointsCount, cache[k + triangleCount]);
                            }
                        }
                    }

                    if(isPrevRightAPWithLeft || (pl.index == anEndPoint.index && isPrevLeftAPWithRight))
                    {
                        const auto triangleArea = std::fabsl(CM::SignedArea(aStartPoint, r, pr));
                        const auto triangleCount = 1 + PointsInTriangle(aStartPoint, r, pr, someBelowPointsCounts, someCollinearPointsCounts);
                        for(size_t i = startIndices[pr.index]; i < rightPointsWithStart.size(); ++i)
                        {
                            const auto& qr = rightPointsWithStart[i];
                            assert(qr.index != r.index && qr.index != pr.index);
                            if (CM::Orientation(qr, pr, r) >= CM::ORIENTATION::COLLINEAR &&
                                CM::Orientation(qr, pr, aStartPoint) >= CM::ORIENTATION::COLLINEAR)
                            {
                                locInsertOrUpdateMinArea<false>(l.index, pl.index, pr.index, qr.index, r.index,
                                                                currArea + triangleArea, aMaxArea, aTotalPointsCount, cache[k + triangleCount]);
                            }
                        }
                    }
                }
            }

            result.myPointsCount += (result.myPointsCount > 0) ? 2 : 0;
            return result;
        }

        void locClearCache(std::vector<std::unordered_map<size_t, Entry>>& anOutCache)
        {
            for(size_t i = 0; i < anOutCache.size(); ++i)
            {
                anOutCache[i].clear();
            }
        }

        void locReconstructHull(const std::vector<std::unordered_map<size_t, Entry>>& aCache,
                                const std::vector<CM::Point2>& somePoints,
                                const std::vector<std::vector<int>>& someBelowPointsCounts,
                                const std::vector<std::vector<int>>& someCollinearPointsCounts,
                                const CM::Point2& aStartPoint,
                                const CM::Point2& anEndPoint,
                                const size_t anOptimalCount,
                                std::vector<size_t>& someOutHullIndices)
        {
            Entry entry;
            entry.a = std::numeric_limits<long double>::infinity();
            size_t k = anOptimalCount;
            for(const auto& [key, value] : aCache[k])
            {
                if(value.pl == anEndPoint.index && value.pr == aStartPoint.index)
                {
                    if(value.a < entry.a) entry = value;
                }
            }

            std::vector<size_t> left = { entry.pl }, right = { entry.pr };
            while(k > 0)
            {
                if(entry.isLeft)
                {
                    k -= 1 + PointsInTriangle(anEndPoint, somePoints[entry.l], somePoints[entry.q], someBelowPointsCounts, someCollinearPointsCounts);
                    left.emplace_back(entry.l);
                    entry = aCache[k].find(KEY(entry.q, entry.l, entry.r, entry.pr, somePoints.size()))->second;
                }
                else
                {
                    k -= 1 + PointsInTriangle(aStartPoint, somePoints[entry.r], somePoints[entry.q], someBelowPointsCounts, someCollinearPointsCounts);
                    right.emplace_back(entry.r);
                    entry = aCache[k].find(KEY(entry.l, entry.pl, entry.q, entry.r, somePoints.size()))->second;
                }
            }

            someOutHullIndices = {};
            someOutHullIndices.resize(left.size() + right.size());
            std::copy(left.begin(), left.end(), someOutHullIndices.begin());
            std::copy(right.begin(), right.end(), someOutHullIndices.begin() + left.size());
        }
    }

    AntipodalResult AntipodalAlgorithm(
            const std::vector<CM::Point2>& somePoints,
            size_t aMaxPointsCount,
            long double aMaxArea,
            long double aMaxDiameter,
            bool aShouldReconstructHull)
    {
        AntipodalResult result;
        if(std::min(somePoints.size(), aMaxPointsCount) < 3)
        {
            return result;
        }

        // Create a 2-dimensional array containing the number of points right below any segment
        std::vector<std::vector<int>> pointsBelowCounts(somePoints.size(), std::vector<int>(somePoints.size(), 0));
        std::vector<std::vector<int>> collinearPointsCounts(somePoints.size(), std::vector<int>(somePoints.size(), 0));

        {
            // Create a 2-dimensional array containing for each point the sequence of clockwise sorted points around it
            std::vector<std::vector<CM::Point2>> clockwiseSortedPoints(somePoints.size(), std::vector<CM::Point2>(somePoints.size() - 1));
            {
                for(size_t i = 0; i < somePoints.size(); ++i)
                {
                    std::copy(somePoints.begin(), somePoints.begin() + i, clockwiseSortedPoints[i].begin());
                    std::copy(somePoints.begin() + i + 1, somePoints.end(), clockwiseSortedPoints[i].begin() + i);
                    SortPointsClockWiseAroundPoint(somePoints[i], clockwiseSortedPoints[i]);
                }
            }
            CountPointsBelowAllSegments(somePoints, clockwiseSortedPoints, pointsBelowCounts, collinearPointsCounts);
        }

        const auto maxAllowedDiameter2 = aMaxDiameter * aMaxDiameter;
        std::pair<size_t, size_t> bestPair;
        std::vector<std::unordered_map<size_t, Entry>> cache {somePoints.size(),
                                                              std::unordered_map<size_t, Entry>{}};

        // Process each distinct pair of points having squared diameter at most equal to maxDiameter2
        for (size_t i = 0; i < somePoints.size(); ++i)
        {
            for (size_t j = i + 1; j < somePoints.size(); ++j)
            {
                AntipodalResult pairResult { std::numeric_limits<long double>::infinity(), 0 };
                if(CM::SquaredDistance(somePoints[i], somePoints[j]) <= maxAllowedDiameter2)
                {
                    std::vector<CM::Point2> leftPoints, rightPoints;
                    locPartitionLeftAndRightPoints(somePoints,somePoints[i],somePoints[j],leftPoints,rightPoints);
                    pairResult = locProcessSegment(somePoints[i], somePoints[j], somePoints, leftPoints, rightPoints, pointsBelowCounts, collinearPointsCounts, somePoints.size(), aMaxPointsCount, aMaxArea, cache);
                    if( pairResult.myPointsCount > result.myPointsCount ||
                        (pairResult.myPointsCount == result.myPointsCount && pairResult.myHullArea < result.myHullArea))
                    {
                        result.myPointsCount = pairResult.myPointsCount;
                        result.myHullArea = pairResult.myHullArea;
                        result.myHasFoundSolution = true;
                        bestPair = std::make_pair(i, j);
                    }
                    locClearCache(cache);
                }
#ifdef DEBUG_ANTIPODAL
                result.results.emplace_back(pairResult.myHullArea, pairResult.myPointsCount);
#endif

#ifdef VERBOSE_ANTIPODAL
                std::cout   << "Done pair (" + std::to_string(i) + ", " + std::to_string(j) + "). Area: " <<
                                pairResult.myHullArea << ", Count: " << pairResult.myPointsCount << std::endl;
#endif
            }
        }

        if(result.myHasFoundSolution && aShouldReconstructHull)
        {
            const auto& s = somePoints[bestPair.first];
            const auto& t = somePoints[bestPair.second];
            std::vector<CM::Point2> leftPoints, rightPoints;
            locPartitionLeftAndRightPoints(somePoints, s, t, leftPoints,rightPoints);
            locProcessSegment(s, t, somePoints, leftPoints, rightPoints, pointsBelowCounts, collinearPointsCounts, somePoints.size(), aMaxPointsCount, aMaxArea, cache);
            locReconstructHull(cache, somePoints,pointsBelowCounts,
                               collinearPointsCounts, s, t,
                               result.myPointsCount - 2, result.myHullIndices);
        }

        return result;
    }
}