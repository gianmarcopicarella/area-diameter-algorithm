//
// Created by Gianmarco Picarella on 27/12/24.
//

#include "Eppstein.h"
#include "CustomMath.h"
#include "Utils.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cassert>

#define MOD(x, n) ((x) % (n))

namespace MT
{
    namespace
    {
        constexpr size_t locIndex(
                const size_t i1, const size_t i2,
                const size_t i3, const size_t i4,
                const std::array<size_t, 4>& axis)
        {
            return i1 * axis[0] + i2 * axis[1] + i3 * axis[2] + i4 * axis[3];
        }

        constexpr int locPointsWithinTriangleCount(
                const std::vector<std::vector<int>>& somePointsBelowCounts,
                const std::vector<std::vector<int>>& somePointsCollinearCounts,
                const int m,
                const CM::Point2& pi,
                const CM::Point2& pj,
                const CM::Point2& pl)
        {
            const auto baseCount =
                    somePointsCollinearCounts[pi.index][pj.index] +
                    somePointsCollinearCounts[pj.index][pl.index] +
                    PointsInTriangle(pi, pj, pl, somePointsBelowCounts, somePointsCollinearCounts);
            if(m == 3 && !CM::AreCollinear(pi, pj, pl))
            {
                return baseCount + somePointsCollinearCounts[pl.index][pi.index];
            }
            return baseCount;
        }

        int locGetFirstClockWiseUpIndex(const std::vector<CM::Point2>& somePoints, const CM::Point2& aRefPoint)
        {
            for(int i = 0; i < somePoints.size() - 1; ++i)
            {
                if(somePoints[i].y < aRefPoint.y && somePoints[i+1].y >= aRefPoint.y)
                {
                    return i + 1;
                }
            }
            return 0;
        }

        void locGetFirstLeftAndRight(
                const std::vector<CM::Point2>& somePoints,
                const CM::Point2& aStartPoint,
                const CM::Point2& anEndPoint,
                std::vector<size_t>& someOutSortedIndicesBySlope)
        {
            using PointWrapper = std::tuple<size_t, CM::ORIENTATION, long double, long double>;
            std::vector<PointWrapper> filteredPoints;
            filteredPoints.reserve(somePoints.size());
            for (size_t i = 0; i < somePoints.size(); ++i)
            {
                if (somePoints[i].y >= aStartPoint.y)
                {
                    const auto angle = CM::Angle(anEndPoint, somePoints[i]);
                    const auto dist = CM::SquaredDistance(anEndPoint, somePoints[i]);
                    const auto orientation = CM::Orientation(aStartPoint, anEndPoint, somePoints[i]);
                    if (orientation == CM::ORIENTATION::COLLINEAR)
                    {
                        if (somePoints[i].y > anEndPoint.y)
                        {
                            filteredPoints.emplace_back(i, CM::ORIENTATION::LEFT_TURN, angle, dist);
                        }
                    }
                    else
                    {
                        filteredPoints.emplace_back(i, orientation, angle, dist);
                    }
                }
            }

            // Find sorted sequence of left and right points
            std::vector<size_t> sortedLeft, sortedRight;

            {
                int startIndex = -1;
                for(int i = 0; i < filteredPoints.size(); ++i)
                {
                    const auto orientation = std::get<1>(filteredPoints[i]);
                    const auto angle = std::get<2>(filteredPoints[i]);
                    const auto distance = std::get<3>(filteredPoints[i]);
                    if(orientation == CM::ORIENTATION::LEFT_TURN &&
                       (startIndex == -1 || angle > std::get<2>(filteredPoints[startIndex]) ||
                        (angle == std::get<2>(filteredPoints[startIndex]) && distance > std::get<3>(filteredPoints[startIndex]))))
                    {
                        startIndex = i;
                    }
                }

                if(startIndex > -1)
                {
                    sortedLeft.emplace_back(std::get<0>(filteredPoints[startIndex]));
                    for(int j = MOD(startIndex + 1, filteredPoints.size()); j != startIndex; j = MOD(j + 1, filteredPoints.size()))
                    {
                        const auto orientation = std::get<1>(filteredPoints[j]);
                        if(orientation == CM::ORIENTATION::LEFT_TURN)
                        {
                            sortedLeft.emplace_back(std::get<0>(filteredPoints[j]));
                        }
                    }
                }
            }

            {
                int startIndex = -1;

                for(int i = 0; i < filteredPoints.size(); ++i)
                {
                    const auto orientation = std::get<1>(filteredPoints[i]);
                    if(orientation == CM::ORIENTATION::RIGHT_TURN)
                    {
                        if(startIndex == -1)
                        {
                            startIndex = i;
                        }
                        else
                        {
                            const auto pointIndex = std::get<0>(filteredPoints[i]);
                            const auto distance = std::get<3>(filteredPoints[i]);
                            auto angle = std::get<2>(filteredPoints[i]);

                            if(somePoints[pointIndex].y > anEndPoint.y)
                            {
                                angle += 2.f * M_PI;
                            }

                            const auto otherPointIndex = std::get<0>(filteredPoints[startIndex]);
                            auto otherAngle = std::get<2>(filteredPoints[startIndex]);
                            if(somePoints[otherPointIndex].y > anEndPoint.y)
                            {
                                otherAngle += 2.f * M_PI;
                            }

                            if( angle > otherAngle ||
                                (angle == otherAngle && distance > std::get<3>(filteredPoints[startIndex])))
                            {
                                startIndex = i;
                            }
                        }
                    }
                }

                if(startIndex > -1)
                {
                    sortedRight.emplace_back(std::get<0>(filteredPoints[startIndex]));
                    for(int j = MOD(startIndex + 1, filteredPoints.size()); j != startIndex; j = MOD(j + 1, filteredPoints.size()))
                    {
                        const auto orientation = std::get<1>(filteredPoints[j]);
                        if(orientation == CM::ORIENTATION::RIGHT_TURN)
                        {
                            sortedRight.emplace_back(std::get<0>(filteredPoints[j]));
                        }
                    }
                }
            }

            std::vector<size_t> test;
            test.reserve(filteredPoints.size());

            std::merge(sortedLeft.begin(), sortedLeft.end(),
                       sortedRight.begin(), sortedRight.end(),
                       std::back_inserter(test), [&](const auto& aRightPointIndex, const auto& aLeftPointIndex){
                const auto& oppositePoint = CM::Point2 {
                        2.f * anEndPoint.x - somePoints[aRightPointIndex].x,
                        2.f * anEndPoint.y - somePoints[aRightPointIndex].y,
                        somePoints[aRightPointIndex].index
                };
                return ArePointsClockwise(anEndPoint, oppositePoint, somePoints[aLeftPointIndex]);
            });

            // Merge left and right
            size_t i = 0, j = 0;
            while(i < sortedLeft.size() && j < sortedRight.size())
            {
                const auto& oppositePoint = CM::Point2 {
                        2.f * anEndPoint.x - somePoints[sortedRight[j]].x,
                        2.f * anEndPoint.y - somePoints[sortedRight[j]].y,
                        somePoints[sortedRight[j]].index
                };
                if(ArePointsClockwise(anEndPoint, somePoints[sortedLeft[i]], oppositePoint))
                {
                    someOutSortedIndicesBySlope.emplace_back(sortedLeft[i++]);
                }
                else
                {
                    someOutSortedIndicesBySlope.emplace_back(sortedRight[j++]);
                }
            }

            if(i >= sortedLeft.size())
            {
                while(j < sortedRight.size())
                {
                    someOutSortedIndicesBySlope.emplace_back(sortedRight[j++]);
                }
            }
            else
            {
                while(i < sortedLeft.size())
                {
                    someOutSortedIndicesBySlope.emplace_back(sortedLeft[i++]);
                }
            }

            assert(test.size() == someOutSortedIndicesBySlope.size());
            for(auto i = 0; i < test.size(); ++i)
            {
                assert(test[i] == someOutSortedIndicesBySlope[i]);
            }
        }
    }

    EppsteinResult EppsteinAlgorithm(
            const std::vector<CM::Point2>& somePoints,
            const int aMaxPointsCount,
            const float aMaxArea,
            const bool aShouldReconstructHull)
    {
        // Create a 2-dimensional array containing for each point the sequence of clockwise sorted points around it
        const size_t pointsCount = somePoints.size();
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
        std::vector<std::vector<int>> pointsBelowCounts(pointsCount, std::vector<int>(pointsCount, 0));
        std::vector<std::vector<int>> collinearPointsCounts(pointsCount, std::vector<int>(pointsCount, 0));
        CountPointsBelowAllSegments(somePoints, clockwiseSortedPoints, pointsBelowCounts, collinearPointsCounts);

        // Sort points by y-coordinate
        std::vector<CM::Point2> sortedPoints(somePoints);
        std::sort(sortedPoints.begin(), sortedPoints.end(), CM::SortPointsVertically);

        // Create a 4-dimensional array storing the minimum areas
        std::vector<long double> minimumAreas((aMaxPointsCount + 1) * pointsCount * pointsCount * pointsCount,
                                        std::numeric_limits<long double>::infinity());
        std::array<size_t, 4> axis{pointsCount * pointsCount * pointsCount, pointsCount * pointsCount, pointsCount, 1};

#ifdef DEBUG_EPPSTEIN
        std::vector<long double> results(aMaxPointsCount + 1, std::numeric_limits<long double>::infinity());
        std::string txt;
#endif

        for (const auto& pi : sortedPoints)
        {
            std::fill_n(&minimumAreas[0] + locIndex(2, pi.index, 0, 0, axis), axis[1], 0);
            for (int m = 3; m < aMaxPointsCount + 1; ++m)
            {
                const auto& clockWisePoints = clockwiseSortedPoints[pi.index];
                const auto startIndex = locGetFirstClockWiseUpIndex(clockWisePoints, pi);
                for (int j = startIndex, count = 0; count < (pointsCount-1) && clockWisePoints[j].y >= pi.y; ++count, j = MOD(j + 1, pointsCount - 1))
                {
                    const auto& pj = clockWisePoints[j];
                    const auto& clockWisePointsAbove = clockwiseSortedPoints[pj.index];
                    std::vector<size_t> sortedIndicesBySlope;
                    locGetFirstLeftAndRight(clockWisePointsAbove, pi, pj, sortedIndicesBySlope);
                    auto minArea = std::numeric_limits<long double>::infinity();

#ifdef DEBUG_EPPSTEIN
                    if(m ==3)
                    {
                        txt += "(" + std::to_string(pi.index) + " | " + std::to_string(pj.index) + "): \n";
                    }
                    auto debug_counter = 0;
#endif
                    for(const auto l : sortedIndicesBySlope)
                    {
                        const auto& pl = clockWisePointsAbove[l];

#ifdef DEBUG_EPPSTEIN
                        if(m == 3)
                        {
                            txt += std::to_string(pl.index);
                            if(debug_counter++ < sortedIndicesBySlope.size()-1)
                            {
                                txt += ", ";
                            }
                            else
                            {
                                txt += "\n";
                            }
                        }
#endif
                        const auto pointsInTriangleCount = 1 + locPointsWithinTriangleCount(pointsBelowCounts, collinearPointsCounts, m, pi, pj, pl);
                        if(pointsInTriangleCount <= m && CM::Orientation(pi, pj, pl) >= CM::ORIENTATION::COLLINEAR)
                        {
                            const auto currentArea =
                                    minimumAreas[locIndex(m - pointsInTriangleCount, pi.index, pl.index, pj.index, axis)] +
                                    std::fabs(CM::SignedArea(pi, pj, pl));
                            if(currentArea < minArea && currentArea <= aMaxArea)
                            {
                                minArea = currentArea;
                            }
                        }
#ifdef DEBUG_EPPSTEIN
                        results[m - 3] = std::min(results[m - 3], minArea);
#endif
                        minimumAreas[locIndex(m, pi.index, pj.index, pl.index, axis)] = minArea;
                    }
                }
            }
        }

        EppsteinResult result;

#ifdef DEBUG_EPPSTEIN
        std::ofstream myfile;
        myfile.open ("slopes.txt");
        myfile << txt;
        myfile.close();
        myfile.open ("results.txt");
        for (const auto r: results)
        {
            myfile << std::fixed << std::setprecision(8) << std::to_string(r) << ", ";
        }
        myfile.close();
        result.results = results;
#endif

        return result;
    }
}