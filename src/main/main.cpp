#include "raylib.h"

#include "../common/Eppstein.h"
#include "../common/Parser.h"

#include <cassert>
#include <iostream>
#include <iomanip>
#include <chrono>

#include <filesystem>
namespace fs = std::filesystem;

using namespace MT;

size_t CountFilesInDirectory(const std::string& aPath)
{
    size_t count = 0;
    for (const auto& entry : fs::directory_iterator(aPath))
    {
        ++count;
    }
    return count;
}

bool AreThereCollinearPoints(const std::vector<CM::Point2>& somePoints)
{
    for(const auto& pi : somePoints)
    {
        for (const auto& pj: somePoints)
        {
            for (const auto& pl: somePoints)
            {
                if (pi.index != pj.index && pi.index != pl.index && pj.index != pl.index)
                {
                    if(CM::AreCollinear(pi, pj, pl))
                    {
                        return true;
                    }
                }
            }
        }
    }
    return false;
}

int main(void)
{
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::seconds;

    const std::string samplesPath = "../../data/samples";
    for (const auto & entry : fs::directory_iterator(samplesPath))
    {
        if(entry.is_directory())
        {
            const auto filesCount = CountFilesInDirectory(entry.path());
            std::cout << "TESTING SAMPLES IN FOLDER: " << entry.path() << std::endl;
            for(size_t i = 0; i < filesCount / 2; ++i)
            {
                const auto sampleFilename = std::string("data_") + std::to_string(i) + ".txt";
                const auto resultsFilename =std::string("results_") + std::to_string(i) + ".txt";
                std::vector<CM::Point2> points;
                std::vector<long double> solutionsPY;
                SZ::ReadPointsFromFile(entry.path() / fs::path{sampleFilename}, points);
                const auto areThereCollinearPoints = AreThereCollinearPoints(points);
                constexpr auto maxAllowedArea = 250;//std::numeric_limits<long double>::infinity();
                constexpr auto shouldReconstructHull = true;

                auto t1 = high_resolution_clock::now();
                const auto res = MT::EppsteinAlgorithm(points, 80, maxAllowedArea, shouldReconstructHull);
                auto t2 = high_resolution_clock::now();

                auto maxSolutionsDistance = 0.l;
/*
                if(SZ::ReadSolutionsFromFile(entry.path() / fs::path{resultsFilename}, solutionsPY))
                {
                    assert(solutionsPY.size() == (res.results.size() - 3));
                    for(int m = 0; m < solutionsPY.size(); ++m)
                    {
                        const auto diff = res.results[m] - solutionsPY[m];
                        maxSolutionsDistance = std::max(maxSolutionsDistance, std::fabs(diff));
                        assert(CM::IsCloseToZero(diff, 0.001));
                    }
                }
*/
                for(const auto idx : res.myHullIndices) std::cout << idx << ", ";
                std::cout << std::endl;
                std::cout << "area=" << res.myHullArea << ", count=" << res.myPointsCount << std::endl;

                std::cout   << std::fixed << std::setprecision(8) << "File " << (1+i) << "/" << (filesCount / 2) << " done in " << duration_cast<seconds>(t2 - t1).count() << " seconds. Collinear points="
                            << (areThereCollinearPoints ? "yes" : "no") << ". Max distance=" << std::to_string(maxSolutionsDistance) << std::endl;
            }
        }
    }

    return 0;
}