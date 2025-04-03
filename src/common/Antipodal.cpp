//
// Created by Gianmarco Picarella on 06/01/25.
//

#include <iostream>
#include <queue>

#include "Antipodal.h"
#include "Utils.h"

#define KEY(a, b, c, d, n) ((a) * (n) * (n) * (n) + (b) * (n) * (n) + (c) * (n) + (d))

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
            const auto maximumDistance2 = CM::Distance2(aStartPoint, anEndPoint);
            const CM::Point2 startToEnd { anEndPoint.myX - aStartPoint.myX, anEndPoint.myY - aStartPoint.myY, INVALID_INDEX };
            const CM::Point2 endToStart { aStartPoint.myX - anEndPoint.myX, aStartPoint.myY - anEndPoint.myY, INVALID_INDEX };
            for (const auto& point : somePoints)
            {
                if( point.myIndex != aStartPoint.myIndex &&
                    point.myIndex != anEndPoint.myIndex &&
                    CM::Distance2(aStartPoint, point) <= maximumDistance2 &&
                    CM::Distance2(anEndPoint, point) <= maximumDistance2)
                {
                    const CM::Point2 startToPoint { point.myX - aStartPoint.myX, point.myY - aStartPoint.myY, INVALID_INDEX };
                    const CM::Point2 endToPoint { point.myX - anEndPoint.myX, point.myY - anEndPoint.myY, INVALID_INDEX };
                    if( CM::Dot(startToEnd, startToPoint) >= 0 &&
                        CM::Dot(endToStart, endToPoint) >= 0)
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
            if(pl.myIndex != someOutPoints[anOutLength - 1].myIndex)
            {
                someOutPoints[anOutLength++] = pl;
            }

            someOutPoints[anOutLength++] = l;

            if(s.myIndex != someOutPoints[anOutLength - 1].myIndex)
            {
                someOutPoints[anOutLength++] = s;
            }
            if(pr.myIndex != someOutPoints[anOutLength - 1].myIndex)
            {
                someOutPoints[anOutLength++] = pr;
            }
            if( r.myIndex != someOutPoints[anOutLength - 1].myIndex &&
                r.myIndex != someOutPoints[0].myIndex)
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
                    return  CM::ProjectedDistance2(aStartPoint, anEndPoint, a) <=
                            CM::ProjectedDistance2(aStartPoint, anEndPoint, b);
                }
                else
                {
                    return  CM::ProjectedDistance2(aStartPoint, anEndPoint, a) >=
                            CM::ProjectedDistance2(aStartPoint, anEndPoint, b);
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
                                  const PointsInTriangleCache& aPointsCountCache,
                                  std::vector<std::unordered_map<size_t, Entry>>& someOutCaches)
        {
            const auto triangleArea = std::fabsl(CM::SignedArea(aFirstPoint, aSecondPoint, aThirdPoint));
            const auto triangleCount = 1 + PointsInTriangle(aFirstPoint, aSecondPoint, aThirdPoint, aPointsCountCache);

            if(aCurrentCount + triangleCount >= someOutCaches.size())
            {
                return;
            }

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
                    // std::cout << aCurrentCount << ", " << triangleCount << ", " << (aCurrentCount + triangleCount) << std::endl;
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

        std::optional<ConvexArea> locProcessSegment(  const std::vector<CM::Point2>& someLeftPoints,
                                                      const std::vector<CM::Point2>& someRightPoints,
                                                      const PointsInTriangleCache& aPointsCountCache,
                                                      const size_t aMaxPointsCount,
                                                      const long double aMaxArea,
                                                      std::vector<std::unordered_map<size_t, Entry>>& someOutCaches)
        {
            const auto maxAllowedPoints = std::min(aMaxPointsCount - 2, someRightPoints.size() + someLeftPoints.size() - 4);
            if(maxAllowedPoints == 0)
            {
                return std::nullopt;
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
            const auto segmentLength2 = CM::Distance2(startPoint, endPoint);
            std::array<CM::Point2, 6> currentPoly;
            size_t currentPolyLength = 0;
            std::optional<ConvexArea> resultOpt;

            for (size_t k = 0; k < maxAllowedPoints + 1; ++k)
            {
                for (const auto& pair : someOutCaches[k])
                {
                    const auto& entry = pair.second;
                    const auto& l = someLeftPoints[entry.l];
                    const auto& r = someRightPoints[entry.r];

                    if(CM::Distance2(l, r) > segmentLength2)
                    {
                        continue;
                    }

                    const auto& pl = someLeftPoints[entry.pl];
                    const auto& pr = someRightPoints[entry.pr];
                    const auto& currArea = entry.area;

                    if(pl.myIndex == endPoint.myIndex && pr.myIndex == startPoint.myIndex)
                    {
                        if(k > 0 && (!resultOpt || k > resultOpt->myPointsCount || (resultOpt->myPointsCount == k && currArea < resultOpt->myHullArea)))
                        {
                            resultOpt = { currArea, k };
                        }
                        continue;
                    }

                    bool isLeftAntipodal = false, isRightAntipodal = false;
                    locKeepUniquePolyPoints(startPoint, endPoint, l, pl, r, pr, currentPoly, currentPolyLength);

                    ForAllAntipodalPairs(currentPoly, [&](const size_t aFirstIndex, const size_t aSecondIndex){
                        const auto& firstPoint = currentPoly[aFirstIndex];
                        const auto& secondPoint = currentPoly[aSecondIndex];
                        isLeftAntipodal |=  (firstPoint.myIndex == pl.myIndex && secondPoint.myIndex == r.myIndex) |
                                            (firstPoint.myIndex == r.myIndex && secondPoint.myIndex == pl.myIndex);

                        isRightAntipodal |= (firstPoint.myIndex == pr.myIndex && secondPoint.myIndex == l.myIndex) |
                                            (firstPoint.myIndex == l.myIndex && secondPoint.myIndex == pr.myIndex);

                    }, currentPolyLength);

                    if(isLeftAntipodal || (pr.myIndex == startPoint.myIndex && isRightAntipodal))
                    {
                        constexpr auto STEP_LEFT = true;
                        locAddNextHullPoint<STEP_LEFT>(endPoint, l, pl, someLeftPoints,
                                                       entry, k, totalPointsCount, aMaxArea,
                                                       aPointsCountCache, someOutCaches);
                    }

                    if(isRightAntipodal || (pl.myIndex == endPoint.myIndex && isLeftAntipodal))
                    {
                        constexpr auto STEP_LEFT = false;
                        locAddNextHullPoint<STEP_LEFT>(startPoint, r, pr, someRightPoints,
                                                       entry, k, totalPointsCount, aMaxArea,
                                                       aPointsCountCache, someOutCaches);
                    }
                }
            }

            if(resultOpt)
            {
                resultOpt->myPointsCount += 2;
            }

            return resultOpt;
        }

        size_t locClearCache(std::vector<std::unordered_map<size_t, Entry>>& someOutCaches)
        {
            size_t entriesCount { 0 };
            for(size_t i = 0; i < someOutCaches.size(); ++i)
            {
                entriesCount += someOutCaches[i].size();
                someOutCaches[i].clear();
            }
            return entriesCount;
        }

        void locReconstructHull(const std::vector<std::unordered_map<size_t, Entry>>& someCaches,
                                const std::vector<CM::Point2>& someLeftPoints,
                                const std::vector<CM::Point2>& someRightPoints,
                                const PointsInTriangleCache& aPointsCountCache,
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
                if( someLeftPoints[value.pl].myIndex == endPoint.myIndex &&
                    someRightPoints[value.pr].myIndex == startPoint.myIndex)
                {
                    if(value.area < entry.area) entry = value;
                }
            }

            std::vector<size_t> left = { someLeftPoints[entry.pl].myIndex }, right = { someRightPoints[entry.pr].myIndex };
            while(k > 0)
            {
                if(entry.isLeftSide)
                {
                    k -= 1 + PointsInTriangle(endPoint, someLeftPoints[entry.l], someLeftPoints[entry.prev], aPointsCountCache);
                    left.emplace_back(someLeftPoints[entry.l].myIndex);
                    const auto key = KEY(entry.prev, entry.l, entry.r, entry.pr, totalPointsCount);
                    entry = someCaches[k].find(key)->second;
                }
                else
                {
                    k -= 1 + PointsInTriangle(startPoint, someRightPoints[entry.r], someRightPoints[entry.prev], aPointsCountCache);
                    right.emplace_back(someRightPoints[entry.r].myIndex);
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
                                          const size_t aFirstIndex,
                                          const size_t aSecondIndex,
                                          std::vector<CM::Point2>& someOutLeftPoints,
                                          std::vector<CM::Point2>& someOutRightPoints)
        {
            const auto& startPoint = somePoints[aFirstIndex];
            const auto& endPoint = somePoints[aSecondIndex];
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

    std::optional<ConvexArea> AntipodalAlgorithm(
            const std::vector<CM::Point2>& somePoints,
            size_t aMaxPointsCount,
            long double aMaxArea,
            long double aMaxDiameter,
            bool aShouldReconstructHull,
            bool aShouldEnableOptimizations)
    {
        std::optional<BenchmarkInfo> benchmarkInfo = std::nullopt;
        return AntipodalAlgorithmWithBenchmarkInfo(somePoints, benchmarkInfo, aMaxPointsCount, aMaxArea, aMaxDiameter, aShouldReconstructHull, aShouldEnableOptimizations);
    }

    std::optional<ConvexArea> AntipodalAlgorithmWithBenchmarkInfo(
            const std::vector<CM::Point2>& somePoints,
            std::optional<BenchmarkInfo>& anOutBenchmarkInfoOpt,
            size_t aMaxPointsCount,
            long double aMaxArea,
            long double aMaxDiameter,
            bool aShouldReconstructHull,
            bool aShouldEnableOptimizations)
    {
        if(std::min(somePoints.size(), aMaxPointsCount) < 3)
        {
            return std::nullopt;
        }

        // Create a 2-dimensional array containing the number of points right below any segment
        PointsInTriangleCache pointsInTriangleCache;

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
            CountPointsBelowAllSegments(somePoints, clockwiseSortedPoints, pointsInTriangleCache);
        }

        std::vector<std::unordered_map<size_t, Entry>> cache {std::min(somePoints.size(), aMaxPointsCount) - 1,
                                                              std::unordered_map<size_t, Entry>{}};
        std::optional<ConvexArea> resultOpt;
        const auto maxAllowedDiameter2 = aMaxDiameter * aMaxDiameter;

        if(anOutBenchmarkInfoOpt)
        {
            anOutBenchmarkInfoOpt->myAllocatedEntriesCount = 0;
            anOutBenchmarkInfoOpt->myRequiredEntriesCount = 0;
        }

        // Process each distinct pair of points having squared diameter at most equal to maxDiameter2
        for (size_t i = 0; i < somePoints.size(); ++i)
        {
            // std::cout << "[" << i << "]" << std::endl;
            for (size_t j = i + 1; j < somePoints.size(); ++j)
            {
                std::optional<ConvexArea> pairResultOpt;
                if(CM::Distance2(somePoints[i], somePoints[j]) <= maxAllowedDiameter2)
                {
                    std::vector<CM::Point2> leftPoints, rightPoints;
                    locPrepareLeftAndRightPoints(somePoints, i, j, leftPoints, rightPoints);

                    if(aShouldEnableOptimizations && resultOpt && (leftPoints.size() + rightPoints.size() - 2) < resultOpt->myPointsCount)
                    {
                        continue;
                    }

                    //std::cout << "(" << j << ")" << std::endl;
                    // std::cout << "count: " << leftPoints.size() + rightPoints.size() - 2 << std::endl;
                    pairResultOpt = locProcessSegment(leftPoints, rightPoints, pointsInTriangleCache, aMaxPointsCount, aMaxArea, cache);
                    if( pairResultOpt && pairResultOpt->myPointsCount <= aMaxPointsCount &&
                        (!resultOpt || pairResultOpt->myPointsCount > resultOpt->myPointsCount ||
                         (pairResultOpt->myPointsCount == resultOpt->myPointsCount && pairResultOpt->myHullArea < resultOpt->myHullArea)))
                    {
                        resultOpt = { pairResultOpt->myHullArea, pairResultOpt->myPointsCount, Diameter{somePoints[i].myIndex, somePoints[j].myIndex} };
                    }

                    const auto cacheEntriesCount = locClearCache(cache);
                    if(anOutBenchmarkInfoOpt)
                    {
                        anOutBenchmarkInfoOpt->myAllocatedEntriesCount += cacheEntriesCount;
                        anOutBenchmarkInfoOpt->myRequiredEntriesCount += cacheEntriesCount;
                    }
                }
            }
        }

        if(resultOpt && aShouldReconstructHull)
        {
            const auto& diameter = *resultOpt->myDiameterOpt;
            std::vector<CM::Point2> leftPoints, rightPoints;
            locPrepareLeftAndRightPoints(somePoints, diameter.myFirstIndex, diameter.mySecondIndex, leftPoints, rightPoints);
            locProcessSegment(leftPoints, rightPoints, pointsInTriangleCache, aMaxPointsCount, aMaxArea, cache);
            locReconstructHull(cache, leftPoints, rightPoints, pointsInTriangleCache, resultOpt->myPointsCount - 2, resultOpt->myHullIndices);
            if(anOutBenchmarkInfoOpt)
            {
                const auto cacheEntriesCount = locClearCache(cache);
                anOutBenchmarkInfoOpt->myAllocatedEntriesCount += cacheEntriesCount;
                anOutBenchmarkInfoOpt->myRequiredEntriesCount += cacheEntriesCount;
            }
        }

        return resultOpt;
    }
}