//
// Created by Gianmarco Picarella on 30/03/25.
//

#include "AntipodalOptimized.h"

#define KEY(a, b, c, n) ((a) * (n) * (n) + (b) * (n) + (c))

namespace MT
{
    namespace
    {
        enum SIDE
        {
            LEFT = 0,
            RIGHT
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
            size_t o, s, ss;
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
                                                      const long double aMaxArea,
                                                      std::vector<std::unordered_map<size_t, Entry>>& someOutCaches)
        {
            const auto maxAllowedPoints = std::min(aMaxPointsCount - 2, someRightPoints.size() + someLeftPoints.size() - 4);
            if(maxAllowedPoints == 0)
            {
                return std::nullopt;
            }

            for(size_t i = 1; i < someLeftPoints.size(); ++i)
            {
                const auto& pl = someLeftPoints[i];
                const auto key = KEY(0, 0, pl.myIndex, aMaxPointsCount);
                someOutCaches[0].emplace(key, Entry {0, 0, i, 0, SIDE::LEFT});
            }
            for(size_t i = 1; i < someRightPoints.size(); ++i)
            {
                const auto& pr = someRightPoints[i];
                const auto key = KEY(0, 0, pr.myIndex, aMaxPointsCount);
                someOutCaches[0].emplace(key, Entry {0, 0, i, 0, SIDE::RIGHT});
            }

            const auto& startPoint = someLeftPoints[0];
            const auto& endPoint = someLeftPoints.back();
            const auto segmentLength2 = CM::Distance2(startPoint, endPoint);

            std::optional<ConvexArea> resultOpt;

            for (size_t k = 0; k < maxAllowedPoints + 1; ++k)
            {
                for (const auto& pair : someOutCaches[k])
                {
                    const auto& entry = pair.second;
                    if( (entry.side == SIDE::RIGHT && someRightPoints[entry.ss].myIndex == startPoint.myIndex) ||
                        (entry.side == SIDE::LEFT && someLeftPoints[entry.ss].myIndex == endPoint.myIndex))
                    {
                        if(k > 0 && (!resultOpt || k > resultOpt->myPointsCount || (resultOpt->myPointsCount == k && entry.area < resultOpt->myHullArea)))
                        {
                            resultOpt = { entry.area, k };
                        }
                        continue;
                    }

                    if(entry.side == SIDE::RIGHT) // right
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

                        if(k + triangleCount >= someOutCaches.size())
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
                                    const auto key = KEY(l.myIndex, pr.myIndex, someRightPoints[i].myIndex, aMaxPointsCount);
                                    auto nextAreaIter = cache.find(key);
                                    if(nextAreaIter == cache.end())
                                    {
                                        cache.emplace(key, Entry{entry.o, entry.ss, i, currentArea, SIDE::RIGHT});
                                    }
                                    else if(currentArea < nextAreaIter->second.area)
                                    {
                                        nextAreaIter->second.area = currentArea;
                                    }
                                }
                            }
                        }

                        const auto refDistance = CM::PointLinePseudoDistance(r, pr, l);

                        // Left points
                        for(size_t i = entry.o + 1; i < someLeftPoints.size(); ++i) {
                            const auto tempDistance = CM::PointLinePseudoDistance(r, pr, someLeftPoints[i]);
                            if (tempDistance <= refDistance && (l.myIndex == startPoint.myIndex ||
                                                                (CM::Orientation(l, endPoint, someLeftPoints[i]) >=
                                                                 CM::ORIENTATION::COLLINEAR &&
                                                                 CM::Orientation(someLeftPoints[i], l, startPoint) >=
                                                                 CM::ORIENTATION::COLLINEAR)))
                            {
                                if (currentArea <= aMaxArea)
                                {
                                    auto &cache = someOutCaches[k + triangleCount];
                                    const auto key = KEY(pr.myIndex, l.myIndex, someLeftPoints[i].myIndex, aMaxPointsCount);
                                    auto nextAreaIter = cache.find(key);
                                    if (nextAreaIter == cache.end())
                                    {
                                        cache.emplace(key, Entry{entry.ss, entry.o, i, currentArea, SIDE::LEFT});
                                    }
                                    else if (currentArea < nextAreaIter->second.area)
                                    {
                                        nextAreaIter->second.area = currentArea;
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

                        if(k + triangleCount >= someOutCaches.size())
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
                                    const auto key = KEY(r.myIndex, pl.myIndex, someLeftPoints[i].myIndex, aMaxPointsCount);
                                    auto nextAreaIter = cache.find(key);
                                    if(nextAreaIter == cache.end())
                                    {
                                        cache.emplace(key, Entry{entry.o, entry.ss, i, currentArea, SIDE::LEFT});
                                    }
                                    else if(currentArea < nextAreaIter->second.area)
                                    {
                                        nextAreaIter->second.area = currentArea;
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
                                    const auto key = KEY(pl.myIndex, r.myIndex, someRightPoints[i].myIndex, aMaxPointsCount);
                                    auto nextAreaIter = cache.find(key);
                                    if (nextAreaIter == cache.end())
                                    {
                                        cache.emplace(key, Entry{entry.ss, entry.o, i, currentArea, SIDE::RIGHT});
                                    }
                                    else if (currentArea < nextAreaIter->second.area)
                                    {
                                        nextAreaIter->second.area = currentArea;
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
        aMaxPointsCount = std::min(somePoints.size(), aMaxPointsCount);
        if(aMaxPointsCount < 3)
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

        std::vector<std::unordered_map<size_t, Entry>> cache {std::min(somePoints.size(), aMaxPointsCount) + 1,
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
            // Add hull reconstruction here
        }

        return resultOpt;
    }
}