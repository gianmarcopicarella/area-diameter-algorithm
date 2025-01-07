//
// Created by Gianmarco Picarella on 06/01/25.
//

#ifndef MASTER_THESIS_ANTIPODAL_H
#define MASTER_THESIS_ANTIPODAL_H

#include <vector>

#define DEBUG_ANTIPODAL

namespace MT
{
    struct AntipodalResult
    {
        bool myHasFoundSolution = false;
        long double myHullArea = 0;
        size_t myPointsCount = 0;
        std::vector<size_t> myHullIndices;
#ifdef DEBUG_ANTIPODAL
        std::vector<long double> results;
#endif
    };
}


#endif //MASTER_THESIS_ANTIPODAL_H
