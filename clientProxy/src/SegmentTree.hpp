#pragma once 
#include <algorithm>
#include <functional>
#include <numeric>
#include <unordered_map>
#include <set>
#include <vector> 

template <typename T>
class SegmentTreeBase
{
public:
    static_assert(!std::is_unsigned<T>::value, "T must be signed");

    virtual void reset(uint32_t bucketNum) = 0;

    virtual void update(uint32_t left, uint32_t right, const T& val) = 0;

    virtual T getTotalWeight() const = 0;

    virtual uint32_t findIdx(T weight) = 0;

    virtual void clear(uint32_t pos) = 0;

    virtual ~SegmentTreeBase() {}
};


template <typename T>
class ConstantSampling : public SegmentTreeBase<T>
{
private:
    static_assert(std::is_integral<T>::value, "T must be integral");
    std::vector<uint32_t> buckets;
    std::unordered_map<uint32_t, uint32_t> lookUpTable;
public:
    virtual void reset(uint32_t bucketNum) override
    {
        lookUpTable.clear();
        buckets.clear();
    }

    virtual void update(uint32_t left, uint32_t right, const T& val) override
    {
        for (uint32_t i = left; i <= right; ++i)
        {
            if (lookUpTable.count(i))
                continue;
            lookUpTable.insert({ i, buckets.size() });
            buckets.push_back(i);
        }
    }

    virtual T getTotalWeight() const override
    {
        return buckets.size();
    }

    virtual uint32_t findIdx(T weight) override
    {
        return buckets[weight];
    }

    virtual void clear(uint32_t pos) override
    {
        if (buckets.empty())
            return;

        auto it = lookUpTable.find(pos);
        if (it == lookUpTable.end())
            return;
        if (buckets.size() > 1)
        {
            lookUpTable[buckets.back()] = it->second;
            std::swap(buckets.back(), buckets[it->second]);
        }
        lookUpTable.erase(it);
        buckets.pop_back();
    }

    virtual ~ConstantSampling() {}
};


template <typename T>
class SegmentTree : public SegmentTreeBase<T>
{
private:
    uint32_t bucketNum;
    uint32_t lastLayer;
    std::vector<T> tree;
    std::function<T(T, T)> fn;

    inline uint32_t leaf(uint32_t k) const
    {
        k += lastLayer - 1;
        k -= (k >= 2 * bucketNum) * bucketNum;
        return k;
    }

    void update(uint32_t pos, const T& val)
    {
        pos = leaf(pos);
        tree[pos] = fn(tree[pos], val);
        while (pos > 1)
        {
            pos >>= 1;
            tree[pos] = tree[pos << 1] + tree[pos << 1 | 1];
        }
    }

public:
    SegmentTree(uint32_t bucketNum,
        const std::function<T(T, T)>& fn) : bucketNum(bucketNum),
        lastLayer(1 << std::__lg(2 * bucketNum - 1)),
        tree(2 * bucketNum),
        fn(fn)
    {
    }

    virtual void reset(uint32_t bucketNum) override
    {
        this->bucketNum = bucketNum;
        lastLayer = 1 << std::__lg(2 * bucketNum - 1);
        tree.resize(2 * bucketNum);
        std::fill(tree.begin(), tree.end(), 0);
    }

    // set weight[left] ... weight[pos]
    virtual void update(uint32_t left, uint32_t right, const T& val) override
    {
        for (uint32_t i = left; i <= right; ++i)
            update(i, val);
    }

    // find the smallest i such that sum(1, i) > weght
    virtual uint32_t findIdx(T weight) override
    {
        uint32_t curIdx = 1;
        while (curIdx < bucketNum)
        {
            if (tree[curIdx << 1] > weight)
                curIdx <<= 1;
            else
            {
                weight -= tree[curIdx << 1];
                curIdx = curIdx << 1 | 1;
            }
        }
        return curIdx >= lastLayer ? curIdx - lastLayer + 1 : curIdx + bucketNum - lastLayer + 1;
    }

    // clear weight[pos] to zero
    virtual void clear(uint32_t pos) override
    {
        pos = leaf(pos);
        tree[pos] = 0;
        while (pos > 1)
        {
            pos >>= 1;
            tree[pos] = tree[pos << 1] + tree[pos << 1 | 1];
        }
    }

    virtual T getTotalWeight() const
    {
        return tree[1];
    }

    virtual ~SegmentTree() {}
};


template <typename T>
class SegmentTreeLazyPropagation : public SegmentTreeBase<T>
{
private:
    uint32_t bucketNum;
    const T flagInit;
    std::vector<T> tree;
    std::vector<T> flags;
    std::function<T(T, T)> flagUpdate;
    std::function<T(uint32_t, T, T)> sumUpdate;

    inline void pushdown(uint32_t curRoot, uint32_t leftLen, uint32_t rightLen)
    {
        flags[curRoot << 1] = flagUpdate(flags[curRoot << 1], flags[curRoot]);
        flags[curRoot << 1 | 1] = flagUpdate(flags[curRoot << 1 | 1], flags[curRoot]);
        tree[curRoot << 1] = sumUpdate(leftLen, tree[curRoot << 1], flags[curRoot]);
        tree[curRoot << 1 | 1] = sumUpdate(rightLen, tree[curRoot << 1 | 1], flags[curRoot]);
        flags[curRoot] = flagInit;
    }

    void update(uint32_t curRoot, uint32_t curLeft, uint32_t curRight,
        uint32_t targetLeft, uint32_t targetRight,
        const T& val)
    {
        if (curLeft > curRight || curLeft > targetRight || curRight < targetLeft)
            return;
        if (curLeft >= targetLeft && curRight <= targetRight)
        {
            if (curLeft == curRight)
            {
                flags[curRoot] = flagInit;
                tree[curRoot] = val == flagInit ? 0 : sumUpdate(1, tree[curRoot], val);
            }
            else
            {
                flags[curRoot] = flagUpdate(flags[curRoot], val);
                tree[curRoot] = sumUpdate(curRight - curLeft + 1, tree[curRoot], val);
            }
            return;
        }
        uint32_t mid = (curLeft + curRight) >> 1;
        if (flags[curRoot] != flagInit)
        {
            pushdown(curRoot, mid - curLeft + 1, curRight - mid);
        }
        update(curRoot << 1, curLeft, mid, targetLeft, targetRight, val);
        update(curRoot << 1 | 1, mid + 1, curRight, targetLeft, targetRight, val);
        tree[curRoot] = tree[curRoot << 1] + tree[curRoot << 1 | 1];
    }

    uint32_t findIdx(uint32_t curRoot, uint32_t left, uint32_t right, T weight)
    {
        if (left == right)
            return left;
        uint32_t mid = (left + right) >> 1;
        if (flags[curRoot] != flagInit)
        {
            pushdown(curRoot, mid - left + 1, right - mid);
        }
        if (tree[curRoot << 1] > weight)
            return findIdx(curRoot << 1, left, mid, weight);
        else
            return findIdx(curRoot << 1 | 1, mid + 1, right, weight - tree[curRoot << 1]);
    }

public:
    SegmentTreeLazyPropagation(uint32_t bucketNum, T flagInit,
        const std::function<T(T, T)>& fnPoint,
        const std::function<T(uint32_t, T, T)>& fnRange) : bucketNum(bucketNum),
        flagInit(flagInit),
        tree(4 * bucketNum + 1),
        flags(4 * bucketNum + 1, flagInit),
        flagUpdate(fnPoint),
        sumUpdate(fnRange)
    {
    }

    virtual void reset(uint32_t bucketNum) override
    {
        this->bucketNum = bucketNum;
        tree.resize(4 * bucketNum + 1);
        flags.resize(4 * bucketNum + 1);
        std::fill(tree.begin(), tree.end(), 0);
        std::fill(flags.begin(), flags.end(), flagInit);
    }

    // set weight[left] ... weight[pos]
    virtual void update(uint32_t left, uint32_t right, const T& val) override
    {
        std::vector<uint32_t> toRemove(right - left + 1);
        std::iota(toRemove.begin(), toRemove.end(), left);

        update(1, 1, bucketNum, left, right, val);
    }

    // find the smallest i such that sum(1, i) > weght
    virtual uint32_t findIdx(T weight) override
    {
        return findIdx(1, 1, bucketNum, weight);
    }

    // clear weight[pos] to zero
    virtual void clear(uint32_t pos) override
    {
        update(1, 1, bucketNum, pos, pos, flagInit);
    }

    virtual T getTotalWeight() const
    {
        return tree[1];
    }

    virtual ~SegmentTreeLazyPropagation() {}
};
