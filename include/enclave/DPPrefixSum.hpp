#pragma once
#include <math.h> 
#include <vector>
#include <enclave/Dist.hpp>

std::vector<uint32_t> DPPrefixSum(const std::vector<uint32_t>& binLoads,
    const double epsilon, const uint32_t seed,
    const uint32_t binCapacity,
    uint32_t& additiveError)
{
    static Laplace laplace(seed);
    const uint32_t binCnt = binLoads.size();
    const double epsilon_ = epsilon / log2(binCnt);
    std::vector<uint32_t> prefixSum(binCnt + 1, 0);
    std::vector<uint32_t> pSums(ceil(log2(binCnt + 1)), 0);
    std::vector<double> noisySums(ceil(log2(binCnt + 1)), 0);
    uint32_t curSum = 0;
    additiveError = 0;
    // REST:
    for (uint32_t t = 1; t <= binCnt; ++t)
    {
        uint32_t i = 0;
        uint32_t t_ = 1;
        while ((t_ & t) == 0)
        {
            t_ <<= 1;
            ++i;
        }
        auto targetPSum = pSums.begin() + i;
        *targetPSum = binLoads[t - 1];
        curSum += binLoads[t - 1];
        auto itP = pSums.begin();
        auto itN = noisySums.begin();
        for (; itP != targetPSum; ++itP, ++itN)
        {
            *targetPSum += *itP;
            *itP = 0;
            *itN = 0;
        }
        noisySums[i] = *targetPSum +
            laplace.sample(1.0 / epsilon);
        double curNoise = 0;
        for (t_ = 1, i = 0; t_ <= t; t_ <<= 1, ++i)
            if (t & t_)
                curNoise += noisySums[i];
        prefixSum[t] = (uint32_t)curNoise < prefixSum[t - 1] ?
            prefixSum[t - 1] : (uint32_t)curNoise;
        prefixSum[t] = prefixSum[t] > curSum + binCapacity ? curSum + binCapacity : prefixSum[t];
        prefixSum[t] = prefixSum[t] + binCapacity < curSum ? curSum - binCapacity : prefixSum[t];
        additiveError = additiveError > abs((int)curSum - (int)prefixSum[t]) ? additiveError : abs((int)curSum - (int)prefixSum[t]);
    }
    return prefixSum;
}