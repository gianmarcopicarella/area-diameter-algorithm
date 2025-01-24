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
        void locPartitionLeftAndRightPoints(const std::vector<CM::Point2>& somePoints,
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

        struct Entry
        {
            size_t l, pl, r, pr, prev;
            long double area { std::numeric_limits<long double>::infinity() };
            bool isLeftSide;
        };

        template <bool isLeftSide=true>
        void locInsertOrUpdateMinArea(const size_t l,
                                      const size_t pl,
                                      const size_t r,
                                      const size_t pr,
                                      const size_t prev,
                                      const long double aCurrentArea,
                                      const long double aMaxArea,
                                      const size_t aTotalPointsCount,
                                      std::unordered_map<size_t, Entry>& anOutcache)
        {
            if(aCurrentArea <= aMaxArea)
            {
                const auto key = KEY(l, pl, r, pr, aTotalPointsCount);
                auto nextAreaIter = anOutcache.find(key);
                if(nextAreaIter == anOutcache.end())
                {
                    anOutcache.emplace(key, Entry{l, pl, r, pr, prev, aCurrentArea, isLeftSide});
                }
                else if(aCurrentArea < nextAreaIter->second.area)
                {
                    nextAreaIter->second.area = aCurrentArea;
                    nextAreaIter->second.prev = prev;
                    nextAreaIter->second.isLeftSide = isLeftSide;
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

        template <bool isLeftSide=true>
        void locAddNextHullPoint( const CM::Point2& aFirstPoint,
                                  const CM::Point2& aSecondPoint,
                                  const CM::Point2& aThirdPoint,
                                  const std::vector<CM::Point2>& somePoints,
                                  const Entry& anEntry,
                                  const size_t aCurrentCount,
                                  const size_t aTotalPointsCount,
                                  const long double aMaxArea,
                                  const std::vector<std::vector<int>>& someBelowPointsCounts,
                                  const std::vector<std::vector<int>>& someCollinearPointsCounts,
                                  std::vector<std::unordered_map<size_t, Entry>>& someOutCaches)
        {
            const auto triangleArea = std::fabsl(CM::SignedArea(aFirstPoint, aSecondPoint, aThirdPoint));
            const auto triangleCount = 1 + PointsInTriangle(aFirstPoint, aSecondPoint, aThirdPoint, someBelowPointsCounts, someCollinearPointsCounts);
            size_t i;
            if constexpr (isLeftSide)
            {
                i = anEntry.pl + 1;
            }
            else
            {
                i = anEntry.pr + 1;
            }
            for(; i < somePoints.size(); ++i)
            {
                if( CM::Orientation(somePoints[i], aThirdPoint, aSecondPoint) >= CM::ORIENTATION::COLLINEAR &&
                    CM::Orientation(somePoints[i], aThirdPoint, aFirstPoint) >= CM::ORIENTATION::COLLINEAR)
                {
                    auto& cache = someOutCaches[aCurrentCount + triangleCount];
                    if constexpr (isLeftSide)
                    {
                        locInsertOrUpdateMinArea<isLeftSide>(   anEntry.pl, i, anEntry.r, anEntry.pr, anEntry.l,
                                                                anEntry.area + triangleArea, aMaxArea, aTotalPointsCount, cache);
                    }
                    else
                    {
                        locInsertOrUpdateMinArea<isLeftSide>(   anEntry.l, anEntry.pl, anEntry.pr, i, anEntry.r,
                                                                anEntry.area + triangleArea, aMaxArea, aTotalPointsCount, cache);
                    }
                }
            }
        }

        AntipodalResult locProcessSegment(  const std::vector<CM::Point2>& someLeftPoints,
                                            const std::vector<CM::Point2>& someRightPoints,
                                            const std::vector<std::vector<int>>& someBelowPointsCounts,
                                            const std::vector<std::vector<int>>& someCollinearPointsCounts,
                                            const size_t aMaxPointsCount,
                                            const long double aMaxArea,
                                            std::vector<std::unordered_map<size_t, Entry>>& someOutCaches)
        {
            AntipodalResult result;

            const auto maxAllowedPoints = std::min(aMaxPointsCount, someRightPoints.size() + someLeftPoints.size() - 4);
            if(maxAllowedPoints == 0)
            {
                return result;
            }

            const auto totalPointsCount = someLeftPoints.size() + someRightPoints.size() - 2;
            for(size_t i = 1; i < someLeftPoints.size(); ++i)
            {
                for(size_t j = 1; j < someRightPoints.size(); ++j)
                {
                    const auto key = KEY(0, i, 0, j, totalPointsCount);
                    someOutCaches[0].emplace(key, Entry {0, i, 0, j, 0, 0});
                }
            }

            const auto& startPoint = someLeftPoints[0];
            const auto& endPoint = someLeftPoints.back();
            const auto segmentLength2 = CM::SquaredDistance(startPoint, endPoint);
            std::array<CM::Point2, 6> currentPoly;
            size_t currentPolyLength = 0;

            for (size_t k = 0; k < maxAllowedPoints + 1; ++k)
            {
                for (const auto& pair : someOutCaches[k])
                {
                    const auto& entry = pair.second;
                    const auto& l = someLeftPoints[entry.l];
                    const auto& r = someRightPoints[entry.r];

                    if(CM::SquaredDistance(l, r) > segmentLength2)
                    {
                        continue;
                    }

                    const auto& pl = someLeftPoints[entry.pl];
                    const auto& pr = someRightPoints[entry.pr];
                    const auto& currArea = entry.area;

                    if(pl.index == endPoint.index && pr.index == startPoint.index)
                    {
                        if(k > 0 && (k > result.myPointsCount || (result.myPointsCount == k && currArea < result.myHullArea)))
                        {
                            result.myHullArea = currArea;
                            result.myPointsCount = k;
                        }
                        continue;
                    }

                    bool isLeftAntipodal = false, isRightAntipodal = false;
                    locKeepUniquePolyPoints(startPoint, endPoint, l, pl, r, pr, currentPoly, currentPolyLength);
                    AreAntipodalPairs(currentPoly, currentPolyLength, pl.index,
                                      r.index, l.index, pr.index,
                                      isLeftAntipodal, isRightAntipodal);

                    if(isLeftAntipodal || (pr.index == startPoint.index && isRightAntipodal))
                    {
                        constexpr auto STEP_LEFT = true;
                        locAddNextHullPoint<STEP_LEFT>(endPoint, l, pl, someLeftPoints,
                                                       entry, k, totalPointsCount, aMaxArea,
                                                       someBelowPointsCounts, someCollinearPointsCounts, someOutCaches);
                    }

                    if(isRightAntipodal || (pl.index == endPoint.index && isLeftAntipodal))
                    {
                        constexpr auto STEP_LEFT = false;
                        locAddNextHullPoint<STEP_LEFT>(startPoint, r, pr, someRightPoints,
                                                       entry, k, totalPointsCount, aMaxArea,
                                                       someBelowPointsCounts, someCollinearPointsCounts, someOutCaches);
                    }
                }
            }

            result.myPointsCount += (result.myPointsCount > 0) ? 2 : 0;
            return result;
        }

        void locClearCache(std::vector<std::unordered_map<size_t, Entry>>& someOutCaches)
        {
            for(size_t i = 0; i < someOutCaches.size(); ++i)
            {
                someOutCaches[i].clear();
            }
        }

        void locReconstructHull(const std::vector<std::unordered_map<size_t, Entry>>& someCaches,
                                const std::vector<CM::Point2>& someLeftPoints,
                                const std::vector<CM::Point2>& someRightPoints,
                                const std::vector<std::vector<int>>& someBelowPointsCounts,
                                const std::vector<std::vector<int>>& someCollinearPointsCounts,
                                const size_t anOptimalCount,
                                std::vector<size_t>& someOutHullIndices)
        {
            const auto& startPoint = someLeftPoints[0];
            const auto& endPoint = someLeftPoints.back();
            const auto totalPointsCount = someLeftPoints.size() + someRightPoints.size() - 2;
            size_t k = anOptimalCount;
            Entry entry;

            for(const auto& [key, value] : someCaches[k])
            {
                if( someLeftPoints[value.pl].index == endPoint.index &&
                    someRightPoints[value.pr].index == startPoint.index)
                {
                    if(value.area < entry.area) entry = value;
                }
            }

            std::vector<size_t> left = { someLeftPoints[entry.pl].index }, right = { someRightPoints[entry.pr].index };
            while(k > 0)
            {
                if(entry.isLeftSide)
                {
                    k -= 1 + PointsInTriangle(endPoint, someLeftPoints[entry.l], someLeftPoints[entry.prev], someBelowPointsCounts, someCollinearPointsCounts);
                    left.emplace_back(someLeftPoints[entry.l].index);
                    const auto key = KEY(entry.prev, entry.l, entry.r, entry.pr, totalPointsCount);
                    entry = someCaches[k].find(key)->second;
                }
                else
                {
                    k -= 1 + PointsInTriangle(startPoint, someRightPoints[entry.r], someRightPoints[entry.prev], someBelowPointsCounts, someCollinearPointsCounts);
                    right.emplace_back(someRightPoints[entry.r].index);
                    const auto key = KEY(entry.l, entry.pl, entry.prev, entry.r, totalPointsCount);
                    entry = someCaches[k].find(key)->second;
                }
            }

            someOutHullIndices = {};
            someOutHullIndices.resize(left.size() + right.size());
            std::copy(left.begin(), left.end(), someOutHullIndices.begin());
            std::copy(right.begin(), right.end(), someOutHullIndices.begin() + left.size());
        }

        void locPrepareLeftAndRightPoints(const std::vector<CM::Point2>& somePoints,
                                          const size_t aStartPointIndex,
                                          const size_t anEndPointIndex,
                                          std::vector<CM::Point2>& someOutLeftPoints,
                                          std::vector<CM::Point2>& someOutRightPoints)
        {
            const auto& startPoint = somePoints[aStartPointIndex];
            const auto& endPoint = somePoints[anEndPointIndex];
            locPartitionLeftAndRightPoints(somePoints, startPoint, endPoint, someOutLeftPoints, someOutRightPoints);
            someOutLeftPoints.emplace_back(startPoint);
            someOutLeftPoints.emplace_back(endPoint);
            someOutRightPoints.emplace_back(startPoint);
            someOutRightPoints.emplace_back(endPoint);
            constexpr auto SORT_LEFT_ASCENDING = true, SORT_RIGHT_ASCENDING = false;
            locSortPointsByProjectionLength<SORT_LEFT_ASCENDING>(startPoint, endPoint, someOutLeftPoints);
            locSortPointsByProjectionLength<SORT_RIGHT_ASCENDING>(startPoint, endPoint, someOutRightPoints);
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
                AntipodalResult pairResult;
                if(CM::SquaredDistance(somePoints[i], somePoints[j]) <= maxAllowedDiameter2)
                {
                    std::vector<CM::Point2> leftPoints, rightPoints;
                    locPrepareLeftAndRightPoints(somePoints, i, j, leftPoints, rightPoints);
#ifdef OPT_USE_OPTIMAL_SOLUTION
                    if((leftPoints.size() + rightPoints.size() - 2) < result.myPointsCount)
                    {
                        continue;
                    }
#endif
                    pairResult = locProcessSegment(leftPoints, rightPoints, pointsBelowCounts, collinearPointsCounts, aMaxPointsCount, aMaxArea, cache);
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
            std::vector<CM::Point2> leftPoints, rightPoints;
            locPrepareLeftAndRightPoints(somePoints, bestPair.first, bestPair.second, leftPoints, rightPoints);
            locProcessSegment(leftPoints, rightPoints, pointsBelowCounts, collinearPointsCounts, aMaxPointsCount, aMaxArea, cache);
            locReconstructHull(cache, leftPoints, rightPoints, pointsBelowCounts, collinearPointsCounts, result.myPointsCount - 2, result.myHullIndices);
        }

        return result;
    }
}