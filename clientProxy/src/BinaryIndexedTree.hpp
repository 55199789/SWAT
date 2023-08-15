#pragma once 
#include "SegmentTree.hpp"
#include <set>
#include <vector>
template <typename T>
class BinaryIndexedTree: public SegmentTreeBase<T>
{
private:
    uint32_t binNum;
    std::vector<uint32_t> roots;
    std::vector<T> ft1; // ft1[i] = sum_{j=1}&{i} weight[j] - weight[j - 1]
    std::vector<T> ft2; // ft2[i] = sum_{j=1}&{i} (weight[j] - weight[j - 1]) * j

    static inline constexpr uint32_t hole(uint32_t k)
    {
        return k + (k >> 10);
    }

    // set weight[pos] += val
    void update(uint32_t pos, const T& val)
    {
        for (uint32_t i = pos; i <= binNum; i += i & (-i))
        {
            ft1[hole(i)] += val;
            ft2[hole(i)] += val * pos;
        }
    }

public:
    BinaryIndexedTree(uint32_t binNum): binNum(binNum),
        ft1(hole(binNum) + 1), ft2(hole(binNum) + 1)
    {
        while (binNum > 0)
        {
            roots.push_back(binNum & -binNum);
            binNum -= roots.back();
        }
        std::reverse(roots.begin(), roots.end());
    }
    virtual void reset(uint32_t binNum) override
    {
        this->binNum = binNum;
        ft1.resize(hole(binNum) + 1);
        ft2.resize(hole(binNum) + 1);
        std::fill(ft1.begin(), ft1.end(), 0);
        std::fill(ft2.begin(), ft2.end(), 0);
        roots.clear();
        while (binNum > 0)
        {
            roots.push_back(binNum & -binNum);
            binNum -= roots.back();
        }
        std::reverse(roots.begin(), roots.end());
    }

    // increment weight[left] ... weight[pos] by val
    virtual void update(uint32_t left, uint32_t right, const T& val) override
    {
        std::vector<uint32_t> toRemove(right - left + 1);
        std::iota(toRemove.begin(), toRemove.end(), left);

        update(left, val);
        update(right + 1, -val);
    }

    // find the smallest i such that sum(1, i) > weght
    virtual uint32_t findIdx(T weight) override
    {
        uint32_t prevIdx = roots.front();
        T partialSum1 = 0;
        T partialSum2 = 0;
        uint32_t i = 1;

        // locate the target tree among all trees
        for (; i < roots.size(); ++i)
        {
            uint32_t curIdx = prevIdx + 1;
            T prevSum = (prevIdx + 1) * (partialSum1 + ft1[hole(prevIdx)]) - (partialSum2 + ft2[hole(prevIdx)]);
            T curSum = prevSum + partialSum1 + ft1[hole(prevIdx)] + (curIdx + 1) * ft1[hole(curIdx)] - ft2[hole(curIdx)];
            if (weight >= prevSum && weight < curSum) // in between two trees
            {
                return curIdx;
            }
            else if (weight < prevSum) // find the target tree
            {
                break;
            }
            else // go to the next tree
            {
                partialSum1 += ft1[hole(prevIdx)];
                partialSum2 += ft2[hole(prevIdx)];
                prevIdx += roots[i];
            }
        }

        // adjust to the left child
        prevIdx -= ((prevIdx & -prevIdx) >> 1);

        // while we have not reached the leaf node
        while ((prevIdx & 1) == 0)
        {
            uint32_t curIdx = prevIdx + 1;
            T prevSum = (prevIdx + 1) * (partialSum1 + ft1[hole(prevIdx)]) - (partialSum2 + ft2[hole(prevIdx)]);
            T curSum = prevSum + partialSum1 + ft1[hole(prevIdx)] + (curIdx + 1) * ft1[hole(curIdx)] - ft2[hole(curIdx)];

            // in between the right most node of the left child and the left most node of the right child
            if (weight >= prevSum && weight < curSum)
            {
                return curIdx;
            }
            else if (weight < prevSum) // go left child
            {
                prevIdx -= ((prevIdx & -prevIdx) >> 1);
            }
            else // go right->left child
            {
                partialSum1 += ft1[hole(prevIdx)];
                partialSum2 += ft2[hole(prevIdx)];
                prevIdx += ((prevIdx & -prevIdx) >> 1);
            }
        }

        // we have reached the leaf node
        T prevSum = (prevIdx + 1) * (partialSum1 + ft1[hole(prevIdx)]) - (partialSum2 + ft2[hole(prevIdx)]);
        if (prevSum > weight)
            return prevIdx;
        else
            return prevIdx + 1;
    }

    // clear weight[pos] to zero
    virtual void clear(uint32_t pos) override
    {
        // 1. compute weight[pos]
        T val = 0;
        uint32_t p1 = pos;
        uint32_t p2 = pos + 1;
        while (p1)
        {
            val += ft1[hole(p1)];
            p1 -= p1 & (-p1);
        }

        update(pos, pos, -val);
    }

    virtual T getTotalWeight() const override
    {
        T sum = 0;
        for (uint32_t i = binNum; i; i -= i & -i)
        {
            sum += (binNum + 1) * ft1[hole(i)] - ft2[hole(i)];
        }
        return sum;
    }

    virtual ~BinaryIndexedTree() {}
};
