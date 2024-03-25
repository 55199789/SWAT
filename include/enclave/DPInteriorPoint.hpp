#pragma once 

#include <math.h> 
#include <random>
#include <stdint.h>
#include <vector>
#include <enclave/Dist.hpp>

template<typename Key>
Key DPInteriorPoint(const std::vector<Key>& bin, uint32_t load, double epsilon, uint32_t seed)
{
    static std::mt19937 gen(seed);
    const double baseExp = exp(epsilon);
    std::vector<double> weights(bin.size());
    for (uint32_t i = 0; i < bin.size();i++)
    {
        // in case the bin is full of dummies
        if (i < load)
            weights[i] = fastPower<double>(baseExp, std::min(i, load - i) + 1);
        else
            weights[i] = 0;
    }
    std::discrete_distribution<uint32_t> d(weights.begin(), weights.end());
    return bin[d(gen)];
}