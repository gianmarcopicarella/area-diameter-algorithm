//
// Created by Gianmarco Picarella on 30/03/25.
//

#include "AntipodalOptimized.h"

#define KEY(a, b, c, n) ((a) * (n) * (n) + (b) * (n) + (c))

#include <iostream>

namespace MT
{
    namespace
    {
        enum SIDE
        {
            LEFT_FROM_LEFT = 0,
            LEFT_FROM_RIGHT,
            RIGHT_FROM_RIGHT,
            RIGHT_FROM_LEFT
        };

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

        struct Entry
        {
            size_t o, s, ss, prev;
            long double area { std::numeric_limits<long double>::infinity() };
            SIDE side;
        };

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

        std::optional<ConvexArea> locProcessSegment(  const std::vector<CM::Point2>& someLeftPoints,
                                                      const std::vector<CM::Point2>& someRightPoints,
                                                      const PointsInTriangleCache& aPointsCountCache,
                                                      const size_t aMaxPointsCount,
                                                      const size_t aTotalPointsCount,
                                                      const long double aMaxArea,
                                                      std::vector<std::unordered_map<size_t, Entry>>& someOutCaches)
        {
            const auto maxAllowedPoints = std::min(aMaxPointsCount - 2, someRightPoints.size() + someLeftPoints.size() - 4);
            if(maxAllowedPoints == 0)
            {
                return std::nullopt;
            }

            const auto& startPoint = someLeftPoints[0];
            const auto& endPoint = someLeftPoints.back();
            const auto segmentLength2 = CM::Distance2(startPoint, endPoint);

            std::optional<ConvexArea> resultOpt;

            for(size_t i = 1; i < someLeftPoints.size(); ++i)
            {
                const auto& pl = someLeftPoints[i];
                const auto key = KEY(startPoint.myIndex, endPoint.myIndex, pl.myIndex, aTotalPointsCount);
                someOutCaches[0].emplace(key, Entry {0, 0, i, INVALID_INDEX, 0, SIDE::LEFT_FROM_LEFT});
            }
            for(size_t i = 1; i < someRightPoints.size(); ++i)
            {
                const auto& pr = someRightPoints[i];
                const auto key = KEY(endPoint.myIndex, startPoint.myIndex, pr.myIndex, aTotalPointsCount);
                someOutCaches[0].emplace(key, Entry {0, 0, i, INVALID_INDEX, 0, SIDE::RIGHT_FROM_RIGHT});
            }

            for (size_t k = 0; k < maxAllowedPoints + 1; ++k)
            {
                for (const auto& pair : someOutCaches[k])
                {
                    const auto& entry = pair.second;
                    if( (entry.side >= SIDE::RIGHT_FROM_RIGHT && someRightPoints[entry.ss].myIndex == startPoint.myIndex) ||
                        (entry.side <= SIDE::LEFT_FROM_RIGHT && someLeftPoints[entry.ss].myIndex == endPoint.myIndex))
                    {
                        if(k > 0 && (!resultOpt || k > resultOpt->myPointsCount || (resultOpt->myPointsCount == k && entry.area < resultOpt->myHullArea)))
                        {
                            resultOpt = { entry.area, k };
                        }
                        continue;
                    }

                    if(entry.side >= SIDE::RIGHT_FROM_RIGHT) // right
                    {
                        const auto& l = someLeftPoints[entry.o];
                        const auto& r = someRightPoints[entry.s];
                        const auto& pr = someRightPoints[entry.ss];

                        if(CM::Distance2(l, pr) > segmentLength2)
                        {
                            continue;
                        }

                        const auto triangleArea = std::fabsl(CM::SignedArea(startPoint, r, pr));
                        const auto triangleCount = 1 + PointsInTriangle(startPoint, r, pr, aPointsCountCache);
                        const auto currentArea = entry.area + triangleArea;

                        if(triangleCount + k >= someOutCaches.size())
                        {
                            continue;
                        }

                        // Right points
                        for(size_t i = entry.ss + 1; i < someRightPoints.size(); ++i)
                        {
                            if( CM::Orientation(someRightPoints[i], pr, r) >= CM::ORIENTATION::COLLINEAR &&
                                CM::Orientation(someRightPoints[i], pr, startPoint) >= CM::ORIENTATION::COLLINEAR)
                            {
                                if(currentArea <= aMaxArea)
                                {
                                    auto& cache = someOutCaches[k + triangleCount];
                                    const auto key = KEY(l.myIndex, pr.myIndex, someRightPoints[i].myIndex, aTotalPointsCount);
                                    auto nextAreaIter = cache.find(key);
                                    if(nextAreaIter == cache.end())
                                    {
                                        cache.emplace(key, Entry{entry.o, entry.ss, i, entry.s, currentArea, SIDE::RIGHT_FROM_RIGHT});
                                    }
                                    else if(currentArea < nextAreaIter->second.area)
                                    {
                                        nextAreaIter->second.area = currentArea;
                                        nextAreaIter->second.prev = entry.s;
                                        nextAreaIter->second.side = SIDE::RIGHT_FROM_RIGHT;
                                    }
                                }
                            }
                        }

                        const auto refDistance = CM::PointLinePseudoDistance(r, pr, l);

                        // Left points
                        for(size_t i = entry.o + 1; i < someLeftPoints.size(); ++i) {
                            const auto tempDistance = CM::PointLinePseudoDistance(r, pr, someLeftPoints[i]);
                            if (tempDistance <= refDistance && (l.myIndex == startPoint.myIndex ||
                                                                (CM::Orientation(l, endPoint, someLeftPoints[i]) >= CM::ORIENTATION::COLLINEAR &&
                                                                 CM::Orientation(someLeftPoints[i], l, startPoint) >= CM::ORIENTATION::COLLINEAR)))
                            {
                                if (currentArea <= aMaxArea)
                                {
                                    auto &cache = someOutCaches[k + triangleCount];
                                    const auto key = KEY(pr.myIndex, l.myIndex, someLeftPoints[i].myIndex, aTotalPointsCount);
                                    auto nextAreaIter = cache.find(key);
                                    if (nextAreaIter == cache.end())
                                    {
                                        cache.emplace(key, Entry{entry.ss, entry.o, i, entry.s, currentArea, SIDE::LEFT_FROM_RIGHT});
                                    }
                                    else if (currentArea < nextAreaIter->second.area)
                                    {
                                        nextAreaIter->second.area = currentArea;
                                        nextAreaIter->second.prev = entry.s;
                                        nextAreaIter->second.side = SIDE::LEFT_FROM_RIGHT;
                                    }
                                }
                            }
                        }
                    }
                    else // left
                    {
                        const auto& r = someRightPoints[entry.o];
                        const auto& l = someLeftPoints[entry.s];
                        const auto& pl = someLeftPoints[entry.ss];

                        if(CM::Distance2(r, pl) > segmentLength2)
                        {
                            continue;
                        }

                        const auto triangleArea = std::fabsl(CM::SignedArea(endPoint, l, pl));
                        const auto triangleCount = 1 + PointsInTriangle(endPoint, l, pl, aPointsCountCache);
                        const auto currentArea = entry.area + triangleArea;

                        if(triangleCount + k >= someOutCaches.size())
                        {
                            continue;
                        }

                        // Left points
                        for(size_t i = entry.ss + 1; i < someLeftPoints.size(); ++i)
                        {
                            if( CM::Orientation(someLeftPoints[i], pl, l) >= CM::ORIENTATION::COLLINEAR &&
                                CM::Orientation(someLeftPoints[i], pl, endPoint) >= CM::ORIENTATION::COLLINEAR)
                            {
                                if(currentArea <= aMaxArea)
                                {
                                    auto& cache = someOutCaches[k + triangleCount];
                                    const auto key = KEY(r.myIndex, pl.myIndex, someLeftPoints[i].myIndex, aTotalPointsCount);
                                    auto nextAreaIter = cache.find(key);
                                    if(nextAreaIter == cache.end())
                                    {
                                        cache.emplace(key, Entry{entry.o, entry.ss, i, entry.s, currentArea, SIDE::LEFT_FROM_LEFT});
                                    }
                                    else if(currentArea < nextAreaIter->second.area)
                                    {
                                        nextAreaIter->second.area = currentArea;
                                        nextAreaIter->second.prev = entry.s;
                                        nextAreaIter->second.side = SIDE::LEFT_FROM_LEFT;
                                    }
                                }
                            }
                        }

                        const auto refDistance = CM::PointLinePseudoDistance(l, pl, r);

                        // Right points
                        for(size_t i = entry.o + 1; i < someRightPoints.size(); ++i)
                        {
                            const auto tempDistance = CM::PointLinePseudoDistance(l, pl, someRightPoints[i]);
                            if (tempDistance <= refDistance && (r.myIndex == endPoint.myIndex ||
                                                                (CM::Orientation(r, startPoint, someRightPoints[i]) >=
                                                                 CM::ORIENTATION::COLLINEAR &&
                                                                 CM::Orientation(someRightPoints[i], r, endPoint) >=
                                                                 CM::ORIENTATION::COLLINEAR)))
                            {
                                if (currentArea <= aMaxArea)
                                {
                                    auto &cache = someOutCaches[k + triangleCount];
                                    const auto key = KEY(pl.myIndex, r.myIndex, someRightPoints[i].myIndex, aTotalPointsCount);
                                    auto nextAreaIter = cache.find(key);
                                    if (nextAreaIter == cache.end())
                                    {
                                        cache.emplace(key, Entry{entry.ss, entry.o, i, entry.s, currentArea, SIDE::RIGHT_FROM_LEFT});
                                    }
                                    else if (currentArea < nextAreaIter->second.area)
                                    {
                                        nextAreaIter->second.area = currentArea;
                                        nextAreaIter->second.prev = entry.s;
                                        nextAreaIter->second.side = SIDE::RIGHT_FROM_LEFT;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            if(resultOpt)
            {
                resultOpt->myPointsCount += 2;
            }
            return resultOpt;
        }

        void locReconstructHull(const std::vector<std::unordered_map<size_t, Entry>>& someCaches,
                                const std::vector<CM::Point2>& someLeftPoints,
                                const std::vector<CM::Point2>& someRightPoints,
                                const PointsInTriangleCache& aPointsCountCache,
                                const size_t anOptimalCount,
                                const size_t aTotalPointCount,
                                std::vector<size_t>& someOutHullIndices)
        {
            const auto& startPoint = someLeftPoints[0];
            const auto& endPoint = someLeftPoints.back();
            size_t k = anOptimalCount;
            Entry entry;

            for(const auto& [key, value] : someCaches[k])
            {
                if( (value.side >= SIDE::RIGHT_FROM_RIGHT && someRightPoints[value.ss].myIndex == startPoint.myIndex) ||
                    (value.side <= SIDE::LEFT_FROM_RIGHT && someLeftPoints[value.ss].myIndex == endPoint.myIndex))
                {
                    if(value.area < entry.area) entry = value;
                }
            }

            std::vector<size_t> left = {}, right = {};

            if(entry.side >= SIDE::RIGHT_FROM_RIGHT)
            {
                left.emplace_back(endPoint.myIndex);
            }
            else if(entry.side <= SIDE::LEFT_FROM_RIGHT)
            {
                right.emplace_back(startPoint.myIndex);
            }

            while(k > 0)
            {
                switch (entry.side)
                {
                    case SIDE::LEFT_FROM_LEFT:
                    {
                        left.emplace_back(someLeftPoints[entry.ss].myIndex);
                        k -= 1 + PointsInTriangle(endPoint, someLeftPoints[entry.prev], someLeftPoints[entry.s], aPointsCountCache);
                        if(k == 0)
                        {
                            left.emplace_back(someLeftPoints[entry.s].myIndex);
                            continue;
                        }
                        const auto key = KEY(someRightPoints[entry.o].myIndex, someLeftPoints[entry.prev].myIndex, someLeftPoints[entry.s].myIndex, aTotalPointCount);
                        entry = someCaches[k].find(key)->second;
                        break;
                    }

                    case SIDE::LEFT_FROM_RIGHT:
                    {
                        left.emplace_back(someLeftPoints[entry.ss].myIndex);

                        k -= 1 + PointsInTriangle(startPoint, someRightPoints[entry.prev], someRightPoints[entry.o], aPointsCountCache);
                        if(k == 0)
                        {
                            right.emplace_back(someRightPoints[entry.o].myIndex);
                            continue;
                        }
                        const auto key = KEY(someLeftPoints[entry.s].myIndex, someRightPoints[entry.prev].myIndex, someRightPoints[entry.o].myIndex, aTotalPointCount);
                        entry = someCaches[k].find(key)->second;
                        break;
                    }

                    case SIDE::RIGHT_FROM_RIGHT:
                    {
                        right.emplace_back(someRightPoints[entry.ss].myIndex);
                        k -= 1 + PointsInTriangle(startPoint, someRightPoints[entry.prev], someRightPoints[entry.s], aPointsCountCache);
                        if(k == 0)
                        {
                            right.emplace_back(someRightPoints[entry.s].myIndex);
                            continue;
                        }
                        const auto key = KEY(someLeftPoints[entry.o].myIndex, someRightPoints[entry.prev].myIndex, someRightPoints[entry.s].myIndex, aTotalPointCount);
                        entry = someCaches[k].find(key)->second;
                        break;
                    }

                    case SIDE::RIGHT_FROM_LEFT:
                    {
                        right.emplace_back(someRightPoints[entry.ss].myIndex);

                        k -= 1 + PointsInTriangle(endPoint, someLeftPoints[entry.prev], someLeftPoints[entry.o], aPointsCountCache);
                        if(k == 0)
                        {
                            left.emplace_back(someLeftPoints[entry.o].myIndex);
                            continue;
                        }
                        const auto key = KEY(someRightPoints[entry.s].myIndex, someLeftPoints[entry.prev].myIndex, someLeftPoints[entry.o].myIndex, aTotalPointCount);
                        entry = someCaches[k].find(key)->second;
                        break;
                    }
                }
            }

            someOutHullIndices = {};
            someOutHullIndices.resize(left.size() + right.size());
            std::copy(left.begin(), left.end(), someOutHullIndices.begin());
            std::copy(right.begin(), right.end(), someOutHullIndices.begin() + left.size());
        }
    }


    std::optional<ConvexArea> AntipodalOptimizedAlgorithm(
            const std::vector<CM::Point2>& somePoints,
            size_t aMaxPointsCount,
            long double aMaxArea,
            long double aMaxDiameter,
            bool aShouldReconstructHull,
            bool aShouldEnableOptimizations)
    {
        std::optional<BenchmarkInfo> benchmarkInfo = std::nullopt;
        return AntipodalOptimizedAlgorithmWithBenchmarkInfo(somePoints, benchmarkInfo, aMaxPointsCount, aMaxArea, aMaxDiameter, aShouldReconstructHull, aShouldEnableOptimizations);
    }

    std::optional<ConvexArea> AntipodalOptimizedAlgorithmWithBenchmarkInfo(
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

                    pairResultOpt = locProcessSegment(leftPoints, rightPoints, pointsInTriangleCache, aMaxPointsCount, somePoints.size(), aMaxArea, cache);

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
            locProcessSegment(leftPoints, rightPoints, pointsInTriangleCache, aMaxPointsCount, somePoints.size(), aMaxArea, cache);
            locReconstructHull(cache, leftPoints, rightPoints, pointsInTriangleCache, resultOpt->myPointsCount - 2, somePoints.size(), resultOpt->myHullIndices);
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