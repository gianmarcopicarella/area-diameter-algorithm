//
// Created by Gianmarco Picarella on 20/02/25.
//

#ifndef MASTER_THESIS_TESTINGUTILS_H
#define MASTER_THESIS_TESTINGUTILS_H

enum class Data
{
    SYNTHETIC_UNIFORM = 0,
    SYNTHETIC_GAUSSIAN,
    COUNT
};

enum class Algorithm
{
    EPPSTEIN = 0,
    ANTIPODAL,
    COUNT
};

static benchmark::IterationCount totalAllocatedBytes { 0 };
static bool shouldTrackHeapMemory { false };

static void StartHeapProfiling()
{
    totalAllocatedBytes = 0;
    shouldTrackHeapMemory = true;
}

static void StopHeapProfiling()
{
    shouldTrackHeapMemory = false;
}

void* operator new(size_t sz)
{
    if(shouldTrackHeapMemory)
    {
        totalAllocatedBytes += sz;
    }
    return std::malloc(sz);
}

template <typename T>
T Mean(const std::vector<T>& someValues)
{
    const auto sum = std::accumulate(someValues.begin(), someValues.end(), T { 0 });
    return sum / someValues.size();
}

template <typename T>
T StandardDeviation(const T& aMean, const std::vector<T>& someValues)
{
    T sumOfSquaredResiduals { 0 };
    for(const auto& aValue : someValues)
    {
        sumOfSquaredResiduals += (aMean - aValue) * (aMean - aValue);
    }
    return std::sqrtl(sumOfSquaredResiduals / someValues.size());
}

#endif //MASTER_THESIS_TESTINGUTILS_H
