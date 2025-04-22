#include "../common/Eppstein.h"
#include "../common/Antipodal.h"
#include "../common/AntipodalOptimized.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;


void PointsListToVector(const py::list& somePoints, std::vector<MT::CM::Point2>& someOutPoints)
{
    someOutPoints.reserve(somePoints.size());
    for(size_t i = 0; i < somePoints.size(); ++i)
    {
        const auto& pyPoint = somePoints[i].cast<std::pair<long double, long double>>();
        someOutPoints.emplace_back(pyPoint.first, pyPoint.second, i);
    }
}

py::object AreaOnlyAlgorithmWithCopy(
        const py::list& points,
        const size_t maxPointsCount,
        const long double maxAllowedArea = std::numeric_limits<long double>::infinity(),
        const bool shouldReconstructHull = false,
        const bool shouldEnableOptimizations = false)
{
    std::vector<MT::CM::Point2> wrappedPoints;
    PointsListToVector(points, wrappedPoints);
    const auto& result = MT::EppsteinAlgorithm(wrappedPoints, maxPointsCount, maxAllowedArea, shouldReconstructHull, shouldEnableOptimizations);
    if(result)
    {
        py::tuple pyResult {3};
        pyResult[0] = result->myHullArea;
        pyResult[1] = result->myPointsCount;
        pyResult[2] = py::cast(result->myHullIndices);
        return pyResult;
    }
    else
    {
        return py::none{};
    }
}

template <bool Optimized=false>
py::object AreaDiameterAlgorithmWithCopy(
        const py::list& points,
        const size_t maxPointsCount,
        const long double maxAllowedArea = std::numeric_limits<long double>::infinity(),
        const long double maxAllowedDiameter = std::numeric_limits<long double>::infinity(),
        const bool shouldReconstructHull = false,
        const bool shouldEnableOptimizations = false)
{
    std::vector<MT::CM::Point2> wrappedPoints;
    PointsListToVector(points, wrappedPoints);

    std::optional<MT::ConvexArea> result;
    if constexpr(Optimized)
    {
        result = MT::AntipodalOptimizedAlgorithm(
                wrappedPoints, maxPointsCount, maxAllowedArea, maxAllowedDiameter, shouldReconstructHull, shouldEnableOptimizations);
    }
    else
    {
        result = MT::AntipodalAlgorithm(
                wrappedPoints, maxPointsCount, maxAllowedArea, maxAllowedDiameter, shouldReconstructHull, shouldEnableOptimizations);
    }

    if(result)
    {
        py::tuple pyResult {4};
        pyResult[0] = result->myHullArea;
        pyResult[1] = result->myPointsCount;
        pyResult[2] = py::cast(result->myHullIndices);
        pyResult[3] = py::cast(std::make_tuple(result->myDiameterOpt->myFirstIndex,
                                               result->myDiameterOpt->mySecondIndex));
        return pyResult;
    }
    else
    {
        return py::none{};
    }
}

PYBIND11_MODULE(thesis, module)
{
    module.def("AreaOnly", &AreaOnlyAlgorithmWithCopy, "The Area-only (A) algorithm (Performs an initial copy of the input points)",
          py::arg("points"),
          py::arg("maxPointsCount"),
          py::arg("maxAllowedArea") = std::numeric_limits<long double>::infinity(),
          py::arg("shouldReconstructHull") = false,
          py::arg("shouldEnableOptimizations") = false);

    module.def("AreaDiameter", &AreaDiameterAlgorithmWithCopy<false>, "The Area-Diameter (AD) algorithm (Performs an initial copy of the input points)",
          py::arg("points"),
          py::arg("maxPointsCount"),
          py::arg("maxAllowedArea") = std::numeric_limits<long double>::infinity(),
          py::arg("maxAllowedDiameter") = std::numeric_limits<long double>::infinity(),
          py::arg("shouldReconstructHull") = false,
          py::arg("shouldEnableOptimizations") = false);

    module.def("AreaDiameterOptimized", &AreaDiameterAlgorithmWithCopy<true>, "The optimized Area-Diameter (AD) algorithm (Performs an initial copy of the input points)",
               py::arg("points"),
               py::arg("maxPointsCount"),
               py::arg("maxAllowedArea") = std::numeric_limits<long double>::infinity(),
               py::arg("maxAllowedDiameter") = std::numeric_limits<long double>::infinity(),
               py::arg("shouldReconstructHull") = false,
               py::arg("shouldEnableOptimizations") = false);

    module.attr("__version__") = "dev";
}

