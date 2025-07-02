#include "Parser.h"

namespace MT
{
    void to_json(json& aJson, const Solution& aSolution)
    {
        aJson = json{ {"name", aSolution.myName}, {"max_count", aSolution.myMaxCount} };
        if(aSolution.myMaxArea != std::numeric_limits<long double>::infinity())
        {
            aJson.push_back({"max_area", aSolution.myMaxArea});
        }
        if(aSolution.myMaxDiameter != std::numeric_limits<long double>::infinity())
        {
            aJson.push_back({"max_diameter", aSolution.myMaxDiameter});
        }
        const auto& convexAreaOpt = aSolution.myConvexAreaOpt;
        if(convexAreaOpt)
        {
            aJson["convex_area"] = {
                    { "count", convexAreaOpt->myPointsCount },
                    { "area", convexAreaOpt->myHullArea },
                    { "hull_indices", convexAreaOpt->myHullIndices }
            };
            if(convexAreaOpt->myDiameterOpt)
            {
                std::vector<size_t> indices = {
                        convexAreaOpt->myDiameterOpt->myFirstIndex,
                        convexAreaOpt->myDiameterOpt->mySecondIndex
                };
                aJson["convex_area"].push_back({ "diameter_indices", indices });
            }
        }
    }
    void from_json(const json& aJson, Solution& anOutSolution)
    {
        if(aJson.contains("name"))
        {
            aJson.at("name").get_to(anOutSolution.myName);
        }
        if(aJson.contains("max_area"))
        {
            aJson.at("max_area").get_to(anOutSolution.myMaxArea);
        }
        if(aJson.contains("max_diameter"))
        {
            aJson.at("max_diameter").get_to(anOutSolution.myMaxDiameter);
        }
        if(aJson.contains("max_count"))
        {
            aJson.at("max_count").get_to(anOutSolution.myMaxCount);
        }
        if(aJson.contains("convex_area"))
        {
            const auto& convexArea = aJson.at("convex_area");
            const auto count = convexArea.at("count").get<size_t>();
            const auto area = convexArea.at("area").get<long double>();
            const auto& diameterIndices = convexArea.at("diameter_indices");
            assert(diameterIndices.size() == 2);
            const auto& hullIndices = convexArea.at("hull_indices").get<std::vector<size_t>>();
            assert(hullIndices.size() > 0);
            anOutSolution.myConvexAreaOpt =
                    { area, count, Diameter{ diameterIndices[0].get<size_t>(), diameterIndices[1].get<size_t>() }, hullIndices};
        }
    }

    namespace SZ
    {
        void ReadPointsFromFile(const std::string& aFilePath, std::vector<CM::Point2>& someOutPoints)
        {
            std::ifstream file(aFilePath);
            if(file.is_open())
            {
                const json& data = json::parse(file);
                assert(data.contains("points"));
                assert(data.contains("count"));
                const auto points = data.at("points");
                for(size_t i = 0; i < data.at("count").get<size_t>(); ++i)
                {
                    const auto x = points[i].at("x").get<long double>();
                    const auto y = points[i].at("y").get<long double>();
                    someOutPoints.emplace_back(x, y, i);
                }
            }
        }

        void ReadSolutionsFromFile(const std::string& aFilePath, std::vector<Solution>& someOutSolutions)
        {
            std::ifstream file(aFilePath);
            if(file.is_open())
            {
                const json& data = json::parse(file);
                assert(data.contains("results"));
                assert(data.contains("count"));
                const auto results = data.at("results");
                const auto solutionsCount = data.at("count").get<size_t>();
                someOutSolutions.reserve(solutionsCount);
                for(size_t i = 0; i < solutionsCount; ++i)
                {
                    someOutSolutions.emplace_back(results[i].get<Solution>());
                }
            }
        }
    }
}