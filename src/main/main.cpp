#include "raylib.h"

#include "../common/Eppstein.h"
#include "../common/Parser.h"

#include <cassert>
#include <iostream>
#include <iomanip>

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
                SZ::ReadSolutionsFromFile(entry.path() / fs::path{resultsFilename}, solutionsPY);
                const auto areThereCollinearPoints = AreThereCollinearPoints(points);
                constexpr auto maxAllowedArea = std::numeric_limits<long double>::infinity();
                constexpr auto shouldReconstructHull = true;
                const auto res = MT::EppsteinAlgorithm(points, points.size(), maxAllowedArea, shouldReconstructHull);
                assert(solutionsPY.size() == (res.results.size() - 3));
                auto maxSolutionsDistance = 0.l;
                for(int m = 0; m < solutionsPY.size(); ++m)
                {
                    const auto diff = res.results[m] - solutionsPY[m];
                    maxSolutionsDistance = std::max(maxSolutionsDistance, std::fabs(diff));
                    assert(CM::IsCloseToZero(diff, 0.001));
                }
                std::cout   << std::fixed << std::setprecision(8) << "File " << (1+i) << "/" << (filesCount / 2) << " done. Collinear points="
                            << (areThereCollinearPoints ? "yes" : "no") << ". Max distance=" << std::to_string(maxSolutionsDistance) << std::endl;
            }
        }
    }

    return 0;
}