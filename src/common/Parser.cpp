//
// Created by Gianmarco Picarella on 28/12/24.
//

#include "Parser.h"
#include "../../external/fast-cpp-csv-parser/csv.h"

namespace MT
{
    namespace SZ
    {
        void ReadPointsFromFile(const std::string& aFilePath, std::vector<CM::Point2>& someOutPoints)
        {
            constexpr auto columnsCount = 2;
            io::CSVReader<columnsCount> input(aFilePath);
            input.read_header(io::ignore_extra_column, "x", "y");
            long double x, y;
            size_t index = 0;
            while(input.read_row(x, y))
            {
                someOutPoints.emplace_back(x, y, index++);
            }
        }

        void ReadSolutionsFromFile(const std::string& aFilePath, std::vector<Solution>& someOutSolutions)
        {
            constexpr auto columnsCount = 7;
            io::CSVReader<columnsCount> input(aFilePath);
            input.read_header(io::ignore_extra_column, "max_count", "max_area", "max_diameter", "area", "count", "diameter_point_index_1", "diameter_point_index_2");
            long double ma, md, a;
            size_t mc, c, d1, d2;
            while(input.read_row(mc, ma, md, a, c, d1, d2))
            {
                std::optional<ConvexArea> areaOpt;
                if(c > 0)
                {
                    areaOpt = { a, c, Diameter{ d1, d2 }};
                }
                someOutSolutions.emplace_back(ma, md, mc, areaOpt);
            }
        }

        void WriteLineToCsv(const std::string& aFilePath, const std::function<std::string(void)>& aCsvLineGenerator, bool aShouldAppendToFile)
        {
            std::ofstream file;
            file.open(aFilePath, std::ios::out | (aShouldAppendToFile ? std::ios::app : std::ios::trunc));
            if(file.is_open())
            {
                file << aCsvLineGenerator() << std::endl;
                file.close();
            }
        }
    }
}