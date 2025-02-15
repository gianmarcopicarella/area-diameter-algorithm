//
// Created by Gianmarco Picarella on 27/12/24.
//

#include "Eppstein.h"
#include "CustomMath.h"
#include "Utils.h"

#include <optional>

#define MOD(x, n) ((x) % (n))
#define IDX(m, i, j, l, count) ((m) * (2 * (count) * (count)) + (i) * ((count) * (count)) + (j) * (count) + (l))

// #define EARLY_STOP_OPT

namespace MT
{
    namespace
    {
        constexpr int locPointsWithinTriangleCount(
                const PointsInTriangleCache& aPointsCountCache,
                const size_t m,
                const CM::Point2& pi,
                const CM::Point2& pj,
                const CM::Point2& pl)
        {
            const auto baseCount =
                    aPointsCountCache.myCollinearPointsCount[pi.myIndex][pj.myIndex] +
                    aPointsCountCache.myCollinearPointsCount[pj.myIndex][pl.myIndex] +
                    PointsInTriangle(pi, pj, pl, aPointsCountCache);
            if(m == 3 && !CM::AreCollinear(pi, pj, pl))
            {
                return baseCount + aPointsCountCache.myCollinearPointsCount[pl.myIndex][pi.myIndex];
            }
            return baseCount;
        }

        int locGetFirstClockWiseUpIndex(const std::vector<CM::Point2>& somePoints, const CM::Point2& aRefPoint)
        {
            for(int i = 0; i < somePoints.size() - 1; ++i)
            {
                if(somePoints[i].myY < aRefPoint.myY && somePoints[i+1].myY >= aRefPoint.myY)
                {
                    return i + 1;
                }
            }
            return 0;
        }

        using Wrapper = std::tuple<long double, long double, size_t>;
        void locCollectSortedPoints(
                const std::vector<Wrapper>& somePoints,
                const std::function<bool(const Wrapper&, const Wrapper&)>& aComparator,
                std::vector<size_t>& someOutSortedPoints)
        {
            if(!somePoints.empty())
            {
                const auto iter = std::max_element(somePoints.begin(), somePoints.end(), aComparator);
                const auto iterIndex = static_cast<size_t>(std::distance(somePoints.begin(), iter));
                someOutSortedPoints.emplace_back(std::get<2>(somePoints[iterIndex]));
                for(int j = MOD(iterIndex + 1, somePoints.size()); j != iterIndex; j = MOD(j + 1, somePoints.size()))
                {
                    someOutSortedPoints.emplace_back(std::get<2>(somePoints[j]));
                }
            }
        }

        void locGetFirstLeftAndRight(
                const std::vector<CM::Point2>& somePoints,
                const CM::Point2& aStartPoint,
                const CM::Point2& anEndPoint,
                std::vector<size_t>& someOutSortedIndicesBySlope)
        {
            std::vector<Wrapper> filteredLeft, filteredRight;
            for (size_t i = 0; i < somePoints.size(); ++i)
            {
                if (somePoints[i].myY >= aStartPoint.myY)
                {
                    const auto angle = CM::Angle(anEndPoint, somePoints[i]);
                    const auto dist = CM::Distance2(anEndPoint, somePoints[i]);
                    const auto orientation = CM::Orientation(aStartPoint, anEndPoint, somePoints[i]);
                    if (orientation == CM::ORIENTATION::COLLINEAR)
                    {
                        if (somePoints[i].myY > anEndPoint.myY)
                        {
                            filteredLeft.emplace_back(angle, dist, i);
                        }
                    }
                    else if(orientation == CM::ORIENTATION::RIGHT_TURN)
                    {
                        filteredRight.emplace_back(angle, dist, i);
                    }
                    else
                    {
                        filteredLeft.emplace_back(angle, dist, i);
                    }
                }
            }

            // Find sorted sequence of left and right points
            std::vector<size_t> sortedLeft, sortedRight;
            locCollectSortedPoints(filteredLeft, std::less<>(), sortedLeft);
            locCollectSortedPoints(filteredRight, [&](const auto& a, const auto& b){
                const auto angle = std::get<0>(a) + (somePoints[std::get<2>(a)].myY > anEndPoint.myY ? (2.f * M_PI) : 0);
                const auto otherAngle = std::get<0>(b) + (somePoints[std::get<2>(b)].myY > anEndPoint.myY ? (2.f * M_PI) : 0);
                return angle < otherAngle || (angle == otherAngle && std::get<1>(a) < std::get<1>(b));
            }, sortedRight);

            someOutSortedIndicesBySlope.reserve(sortedLeft.size() + sortedRight.size());
            std::merge(sortedLeft.begin(), sortedLeft.end(),
                       sortedRight.begin(), sortedRight.end(),
                       std::back_inserter(someOutSortedIndicesBySlope), [&](const auto& aRightPointIndex, const auto& aLeftPointIndex){
                        const auto& oppositePoint = CM::Point2 {
                                2.f * anEndPoint.myX - somePoints[aRightPointIndex].myX,
                                2.f * anEndPoint.myY - somePoints[aRightPointIndex].myY,
                                somePoints[aRightPointIndex].myIndex
                        };
                        return ArePointsClockwise(anEndPoint, oppositePoint, somePoints[aLeftPointIndex]);
                    });
        }
    }

    std::optional<ConvexArea> EppsteinAlgorithm(
            const std::vector<CM::Point2>& somePoints,
            const size_t aMaxPointsCount,
            const long double aMaxArea,
            const bool aShouldReconstructHull)
    {
        std::optional<BenchmarkInfo> benchmarkInfo = std::nullopt;
        return EppsteinAlgorithmWithBenchmarkInfo(somePoints, benchmarkInfo, aMaxPointsCount, aMaxArea, aShouldReconstructHull);
    }

    std::optional<ConvexArea> EppsteinAlgorithmWithBenchmarkInfo(
            const std::vector<CM::Point2>& somePoints,
            std::optional<BenchmarkInfo>& anOutBenchmarkInfoOpt,
            size_t aMaxPointsCount,
            long double aMaxArea,
            bool aShouldReconstructHull)
    {
        const auto pointsCount = somePoints.size();
        if(pointsCount < 3 || aMaxPointsCount < 3)
        {
            return std::nullopt;
        }

        // Create a 2-dimensional array containing for each point the sequence of clockwise sorted points around it
        std::vector<std::vector<CM::Point2>> clockwiseSortedPoints(pointsCount, std::vector<CM::Point2>(pointsCount - 1));
        {
            for(int i = 0; i < pointsCount; ++i)
            {
                std::copy(somePoints.begin(), somePoints.begin() + i, clockwiseSortedPoints[i].begin());
                std::copy(somePoints.begin() + i + 1, somePoints.end(), clockwiseSortedPoints[i].begin() + i);
                SortPointsClockWiseAroundPoint(somePoints[i], clockwiseSortedPoints[i]);
            }
        }

        // Create a 2-dimensional array containing the number of points right below any segment
        PointsInTriangleCache pointsCountCache;
        CountPointsBelowAllSegments(somePoints, clockwiseSortedPoints, pointsCountCache);

        // Sort points by y-coordinate
        std::vector<CM::Point2> sortedPoints(somePoints);
        std::sort(sortedPoints.begin(), sortedPoints.end(), CM::SortPointsVertically);

        // Create a 3-dimensional array storing the minimum areas
        const auto maxPointsCount = std::min(aMaxPointsCount, somePoints.size());
        std::vector<long double> minimumAreas((maxPointsCount + 1) * 2 * pointsCount * pointsCount,
                                              std::numeric_limits<long double>::infinity());

        std::optional<ConvexArea> resultOpt;
        std::pair<size_t, size_t> bestIndex;
        size_t minimumAreaIndex { 0 };
        bool hasFoundNewBestIndex { false };

        if(anOutBenchmarkInfoOpt)
        {
            anOutBenchmarkInfoOpt->myCreatedEntriesCount = minimumAreas.size();
            anOutBenchmarkInfoOpt->myMinEntriesCount = 0;
        }

        for (size_t i = 0; i < sortedPoints.size(); ++i)
        {
            const auto& pi = sortedPoints[i];
            if(hasFoundNewBestIndex)
            {
                minimumAreaIndex ^= 1;
                hasFoundNewBestIndex = false;
            }
#ifdef EARLY_STOP_OPT
            if(resultOpt && resultOpt->myPointsCount > (sortedPoints.size() - i))
            {
                break;
            }
#endif
            std::fill_n(&minimumAreas[0] + IDX(2, minimumAreaIndex, 0, 0, pointsCount), pointsCount * pointsCount, 0);
            for (size_t m = 3; m < maxPointsCount + 1; ++m)
            {
                const auto& clockWisePoints = clockwiseSortedPoints[pi.myIndex];
                const auto startIndex = locGetFirstClockWiseUpIndex(clockWisePoints, pi);
                for (int j = startIndex, count = 0; count < (pointsCount-1) && clockWisePoints[j].myY >= pi.myY; ++count, j = MOD(j + 1, pointsCount - 1))
                {
                    const auto& pj = clockWisePoints[j];
                    const auto& clockWisePointsAbove = clockwiseSortedPoints[pj.myIndex];
                    std::vector<size_t> sortedIndicesBySlope;
                    locGetFirstLeftAndRight(clockWisePointsAbove, pi, pj, sortedIndicesBySlope);
                    auto minArea = std::numeric_limits<long double>::infinity();
                    auto lastMinArea = -minArea;

                    for(const auto l : sortedIndicesBySlope)
                    {
                        const auto& pl = clockWisePointsAbove[l];

                        const auto pointsInTriangleCount = 1 + locPointsWithinTriangleCount(pointsCountCache, m, pi, pj, pl);
                        if(pointsInTriangleCount <= m && CM::Orientation(pi, pj, pl) >= CM::ORIENTATION::COLLINEAR)
                        {
                            const auto currentArea =
                                    minimumAreas[IDX(m - pointsInTriangleCount, minimumAreaIndex, pl.myIndex, pj.myIndex, pointsCount)] +
                                    std::fabs(CM::SignedArea(pi, pj, pl));

                            if(currentArea < minArea && currentArea <= aMaxArea)
                            {
                                minArea = currentArea;
                                if(!resultOpt || m > resultOpt->myPointsCount || (m == resultOpt->myPointsCount && minArea < resultOpt->myHullArea))
                                {
                                    hasFoundNewBestIndex |= !resultOpt || bestIndex.first != pi.myIndex;
                                    resultOpt = {minArea, m};
                                    bestIndex = {pi.myIndex, pj.myIndex};
                                }
                            }
                        }

                        if(minArea != lastMinArea && anOutBenchmarkInfoOpt)
                        {
                            lastMinArea = minArea;
                            anOutBenchmarkInfoOpt->myMinEntriesCount += 1;
                        }

                        minimumAreas[IDX(m, minimumAreaIndex, pj.myIndex, pl.myIndex, pointsCount)] = minArea;
                    }
                }
            }
        }

        if(resultOpt && aShouldReconstructHull)
        {
            minimumAreaIndex = hasFoundNewBestIndex ? minimumAreaIndex : (minimumAreaIndex ^ 1);
            size_t m = resultOpt->myPointsCount, j = bestIndex.second;
            resultOpt->myHullIndices.emplace_back(bestIndex.first);
            while (m > 2)
            {
                size_t l { INVALID_INDEX };
                const auto& latestHullPoint = resultOpt->myHullIndices.back();
                const auto& clockWisePointsAbove = clockwiseSortedPoints[j];
                std::vector<size_t> sortedIndicesBySlope;
                locGetFirstLeftAndRight(clockWisePointsAbove, somePoints[bestIndex.first], somePoints[j], sortedIndicesBySlope);
                for (size_t k = sortedIndicesBySlope.size() - 1; k >= 0 ; --k)
                {
                    const auto& current = clockWisePointsAbove[sortedIndicesBySlope[k]];
                    const auto& previous = clockWisePointsAbove[sortedIndicesBySlope[k - 1]];
                    if( (k == 0 ||
                         minimumAreas[IDX(m, minimumAreaIndex, j, current.myIndex, pointsCount)] !=
                         minimumAreas[IDX(m, minimumAreaIndex, j, previous.myIndex, pointsCount)]) &&
                        CM::Orientation(somePoints[latestHullPoint], somePoints[j], current) >= CM::ORIENTATION::COLLINEAR)
                    {
                        l = clockWisePointsAbove[sortedIndicesBySlope[k]].myIndex;
                        break;
                    }
                }
                m -= 1 + locPointsWithinTriangleCount(pointsCountCache, m, somePoints[bestIndex.first], somePoints[j], somePoints[l]);
                resultOpt->myHullIndices.emplace_back(j);
                j = l;
            }
            resultOpt->myHullIndices.emplace_back(j);

            std::vector<CM::Point2> hullPoints;
            GetHullPoints(resultOpt->myHullIndices, somePoints, hullPoints);
            resultOpt->myDiameterOpt = ComputeDiameter(hullPoints);
        }

        return resultOpt;
    }
}