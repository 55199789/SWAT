#pragma once
#include <assert.h>
#include <inttypes.h>
#include <iterator>
#include <list>
#include <memory>
#include <random>
#include <vector>
#include "common/defs.h"
#include "SegmentTree.hpp"
#include "BinaryIndexedTree.hpp" 

template <typename Key>
class QueryEntity
{
public:
    const Key left;
    const Key right;
    const std::chrono::high_resolution_clock::time_point st;
    uint32_t cnt; // number of pending buckets in this query
    bool empty;

    QueryEntity(Key left, Key right)
        : left(left),
        right(right),
        st(std::chrono::high_resolution_clock::now()),
        cnt(0),
        empty(true)
    {
    }

    void notify()
    {
        printf("[%d, %d] completed within %" PRId64 " microseconds\n",
            left, right,
            std::chrono::duration_cast<std::chrono::microseconds>(
                std::chrono::high_resolution_clock::now() - st
            ).count());
    }

    ~QueryEntity()
    {
        // prevent empty queries from being printed
        if (!empty)
        {
            // notify();
        }
    }
};
template <typename Key,
    typename WeightType = int>
class SamplePool
{
private:
    SamplePool() {}

    std::mt19937 gen;

    uint32_t levels;
    uint32_t batchSize;
    double storageOverhead;
    double alpha;
    double delta;
    Key min;
    Key max;
    Key rangeLen;
    std::vector<std::unique_ptr<SegmentTreeBase<WeightType>>> segmentTrees;
    std::vector<std::vector<BucketTag<Key>>> bucketTagLists;
    /**
     *  first layer: level
     *  second layer: buckets in the current level
     *  third layer: pending queries (IDs) in the current bucket
     **/
    std::vector<std::vector<std::vector<typename std::list<QueryEntity<Key>>::iterator>>> invertedIndex;

    std::list<QueryEntity<Key>> queries;

    std::vector<std::vector<std::string>> bucketIDLists;
    std::vector<std::vector<uint32_t>> replicaCntLists;

    std::vector<std::discrete_distribution<uint32_t>> realDistributionLists;
    std::vector<std::discrete_distribution<uint32_t>> fakeDistributionLists;

    std::uniform_real_distribution<double> deltaDist;

public:
    SamplePool(const SamplePool<Key, WeightType>&) = delete;
    SamplePool(SamplePool<Key, WeightType>&&) = delete;
    SamplePool<Key, WeightType>& operator=(const SamplePool<Key, WeightType>&) = delete;
    SamplePool<Key, WeightType>& operator=(SamplePool<Key, WeightType>&&) = delete;

    static SamplePool<Key, WeightType>& getInstance()
    {
        static SamplePool<Key, WeightType> instance;
        return instance;
    }

    void setup(uint32_t levels, uint32_t batchSize, uint32_t seed,
        double storageOverhead, Key min, Key max, Key rangeLen,
        WeightUpdatePolicy policy = WeightUpdatePolicy::Increment)
    {
        gen.seed(seed);
        makeBucketID(0, 0, 0, seed);
        this->levels = levels;
        this->batchSize = batchSize;
        this->storageOverhead = storageOverhead;
        this->min = min;
        this->max = max;
        this->rangeLen = rangeLen;
        alpha = 1.0 / (storageOverhead - 1.0);
        delta = (storageOverhead - 1.0) / storageOverhead;

        invertedIndex.resize(levels);
        bucketTagLists.resize(levels);
        bucketIDLists.resize(levels);
        replicaCntLists.resize(levels);
        realDistributionLists.resize(levels);
        fakeDistributionLists.resize(levels);
        deltaDist = std::uniform_real_distribution<double>(0.0, 1.0);

        for (uint32_t i = 0; i < levels; ++i)
        {
            switch (policy)
            {
            case WeightUpdatePolicy::Increment:
                segmentTrees.emplace_back(
                    new BinaryIndexedTree<WeightType>(0));
                break;
            case WeightUpdatePolicy::Constant:
                segmentTrees.emplace_back(
                    new ConstantSampling<WeightType>()
                );
                break;
            case WeightUpdatePolicy::Multiplication:
                segmentTrees.emplace_back(
                    new SegmentTreeLazyPropagation<WeightType>(
                        0,
                        1,
                        [](int x, int val)
                        { return x * val; },
                        [](uint32_t len, int x, int val)
                        { return x * val; })
                );
                break;
            default:
                throw std::runtime_error("Invalid weight update policy");
            }
        }
    }

    void insertQuery(Key qLeft, Key qRight)
    {
        queries.emplace_back(qLeft, qRight);
        QueryEntity<Key>& query = queries.back();
        auto queryIter = std::prev(queries.end());
        // update weights in sample pool
        bool flag = false;
        for (uint32_t i = 0; i < levels; ++i)
        {
            const auto& bucketTagList = bucketTagLists[i];
            // current layer is empty
            if (bucketTagList.empty())
            {
                continue;
            }
            query.empty = false;

            // find left bound of bucketTagList
            auto leftTag = std::upper_bound(bucketTagList.rbegin(), bucketTagList.rend(), query.left,
                [](const Key& key, const BucketTag<Key>& bucketTag)
                {
                    return bucketTag.right < key;
                });
            if (leftTag == bucketTagList.rbegin())
            {
                continue;
            }

            // find right bound of bucketTagList
            auto rightTag = std::upper_bound(bucketTagList.begin(), bucketTagList.end(), query.right,
                [](const Key& key, const BucketTag<Key>& bucketTag)
                {
                    return bucketTag.left > key;
                });

            // no candidate buckets as all buckets locate at the right of the query
            if (rightTag == bucketTagList.begin())
            {
                continue;
            }

            uint32_t left = bucketTagList.rend() - leftTag + 1;
            uint32_t right = rightTag - bucketTagList.begin();

            // assert(left <= right); // left bound should be smaller than right bound
            if (left > right)
                continue;

            flag = true;
            /**
             * update weights by one of the following approaches:
             *  1. add by 2;
             *  2. set to constant 2;
             *  3. multiply by 2
             * */
            segmentTrees[i]->update(left, right, 2);

            // update the number of pending buckets in this query
            query.cnt += right - left + 1;
            for (uint32_t j = left; j <= right; ++j)
                invertedIndex[i][j - 1].push_back(queryIter);
        }
        if (!flag)
            queries.pop_back();
        // else
        //     printf("insert query:[%d, %d] = %d\n", query.left, query.right, query.cnt);
    }

    void insertQueries(const std::vector<QueryEntity<Key>>& queries)
    {
        for (const auto& query : queries)
        {
            insertQuery(query);
        }
    }

    void sample(std::vector<uint8_t>& buffer)
    {
        if (queries.empty())
            return;
        for (uint32_t i = 0; i < levels; ++i)
        {
            // if current level is empty, skip
            if (bucketTagLists[i].empty())
            {
                continue;
            }

            for (uint32_t j = 0; j < batchSize; ++j)
            {
                // true --> real query, false --> fake query
                const bool queryType = deltaDist(gen) < delta;
                uint32_t id = 0;
                uint32_t bucket = 0;
                if (queryType) // generate real queries
                {
                    WeightType sum = segmentTrees[i]->getTotalWeight();
                    if (sum == 0) // no pending buckets
                    {
                        id = realDistributionLists[i](gen);
                    }

                    else // sample a pending bucket and pick a random replica from it
                    {
                        std::uniform_int_distribution<uint32_t> randWeightGen(0, sum - 1);
                        bucket = segmentTrees[i]->findIdx(randWeightGen(gen)) - 1;
                        segmentTrees[i]->clear(bucket + 1);
                        for (auto& queryIter : invertedIndex[i][bucket])
                        {
                            --queryIter->cnt;
                            if (queryIter->cnt == 0)
                            {
                                queryIter->notify();
                                queries.erase(queryIter);
                            }
                        }
                        invertedIndex[i][bucket].clear();
                        uint32_t prevAccCnt = bucket == 0 ? 0 : replicaCntLists[i][bucket - 1];
                        std::uniform_int_distribution<uint32_t> randBucketGen(1,
                            replicaCntLists[i][bucket] - prevAccCnt);
                        id = randBucketGen(gen) + prevAccCnt - 1;
                    }
                }
                else // generate fake queries
                {
                    id = fakeDistributionLists[i](gen);
                    bucket = std::upper_bound(replicaCntLists[i].begin(),
                        replicaCntLists[i].end(), id) - replicaCntLists[i].begin();
                    if (bucket < bucketTagLists[i].size())
                    {
                        segmentTrees[i]->clear(bucket + 1);
                        for (auto& queryIter : invertedIndex[i][bucket])
                        {
                            --queryIter->cnt;
                            if (queryIter->cnt == 0)
                            {
                                queryIter->notify();
                                queries.erase(queryIter);
                            }
                        }
                        invertedIndex[i][bucket].clear();
                    }
                }
                buffer.insert(buffer.end(), bucketIDLists[i][id].begin(), bucketIDLists[i][id].end());
                buffer.push_back(' ');
            }
        }
    }

    std::vector<uint32_t> update(uint32_t level, std::vector<BucketTag<Key>>& bucketTagList)
    {
        uint32_t bucketCnt = bucketTagList.size();
        for (uint32_t i = 0;i < bucketCnt; ++i)
            if (bucketTagList[i].left > max)
            {
                assert(i > 0);
                bucketCnt = i;
                bucketTagList.resize(bucketCnt);
                break;
            }

        std::vector<std::vector<typename std::list<QueryEntity<Key>>::iterator>>
            curLevelInvertedIndex = invertedIndex[level];
        invertedIndex[level].clear();
        invertedIndex[level].resize(bucketCnt);
        segmentTrees[level]->reset(bucketCnt);

        // transform unqueried buckets in merged levels to the current level
        for (uint32_t i = 0; i <= level; ++i)
        {
            std::vector<std::vector<typename std::list<QueryEntity<Key>>::iterator>>&
                curInvertedIndex = (i == level) ? curLevelInvertedIndex : invertedIndex[i];
            if (curInvertedIndex.empty())
                continue;
            for (uint32_t bucket = 0; bucket < bucketTagLists[i].size(); ++bucket)
            {
                const BucketTag<Key>& curBucketTag = bucketTagLists[i][bucket];
                for (auto& qIter : curInvertedIndex[bucket])
                {
                    Key curLeft = std::max(qIter->left, curBucketTag.left);
                    Key curRight = std::min(qIter->right, curBucketTag.right);


                    // find left bound of bucketTagList
                    auto leftTag = std::upper_bound(bucketTagList.rbegin(), bucketTagList.rend(), curLeft,
                        [](const Key& key, const BucketTag<Key>& bucketTag)
                        {
                            return bucketTag.right < key;
                        });
                    // if (leftTag == bucketTagList.rbegin())
                    // {
                    //     continue;
                    // }

                    uint32_t left = bucketTagList.rend() - leftTag;
                    // find right bound of bucketTagList
                    auto rightTag = std::upper_bound(bucketTagList.begin() + left, bucketTagList.end(), curRight,
                        [](const Key& key, const BucketTag<Key>& bucketTag)
                        {
                            return bucketTag.left > key;
                        });

                    // no candidate buckets as all buckets locate at the right of the query
                    // if (rightTag == bucketTagList.begin() + left)
                    // {
                    //     continue;
                    // }
                    left++;
                    uint32_t right = rightTag - bucketTagList.begin();

                    if (left > right)
                    {
                        qIter->cnt--;
                        if (qIter->cnt == 0)
                        {
                            qIter->notify();
                            queries.erase(qIter);
                        }
                        continue;
                    }
                    assert(left <= right); // left bound should be smaller than right bound

                    /**
                     * update weights by one of the following approaches:
                     *  1. add by 2;
                     *  2. set to constant 2;
                     *  3. multiply by 2
                     * */
                    segmentTrees[level]->update(left, right, 2);

                    for (uint32_t curBucket = left; curBucket <= right; ++curBucket)
                        invertedIndex[level][curBucket - 1].push_back(qIter);
                    qIter->cnt += right - left;
                }
            }
            curInvertedIndex.clear();
            replicaCntLists[i].clear();
            bucketTagLists[i].clear();
            bucketIDLists[i].clear();
        }

        const uint32_t totalBucketCnt = ceil(bucketCnt * storageOverhead);

        std::vector<double> weights(bucketCnt);
        std::vector<double> realWeights(totalBucketCnt);
        std::vector<double> fakeWeights(totalBucketCnt);
        std::vector<uint32_t>& replicaCnt = replicaCntLists[level];
        std::vector<std::string>& bucketIDList = bucketIDLists[level];

        replicaCnt.resize(bucketCnt + 1);
        bucketIDList.resize(totalBucketCnt);
        const Key n = max - min + 1;
        double wSum = 0;
        for (uint32_t i = 0; i < bucketCnt; ++i)
        {
            Key x = std::min(rangeLen + 1, bucketTagList[i].left - min);
            Key y = std::min(rangeLen + 1, max - bucketTagList[i].right);
            weights[i] = 1.0 - \
                (2.0 * bucketTagList[i].left - 1.0 - x) / (2.0 * n - rangeLen) * \
                (1.0 * x / (rangeLen + 1.0)) - \
                (2.0 * n - 2.0 * bucketTagList[i].right + 1.0 - y) / (2.0 * n - rangeLen) * \
                (1.0 * y / (rangeLen + 1.0));
            wSum += weights[i];
        }
        for (auto& w : weights)
            w /= wSum;
        bucketTagLists[level] = std::move(bucketTagList);

        uint32_t curCnt = 0;
        for (uint32_t i = 0; i < bucketCnt; ++i)
        {
            replicaCnt[i] = std::ceil(weights[i] * bucketCnt / alpha);

            std::fill(fakeWeights.begin() + curCnt,
                fakeWeights.begin() + curCnt + replicaCnt[i],
                (alpha / bucketCnt - weights[i] / replicaCnt[i]) / (1.0 / delta - 1.0));
            std::fill(realWeights.begin() + curCnt,
                realWeights.begin() + curCnt + replicaCnt[i],
                weights[i]);

            for (uint32_t j = 0; j < replicaCnt[i]; ++j)
                bucketIDList[curCnt + j] = makeBucketID(level, i, j);
            curCnt += replicaCnt[i];
        }
        replicaCnt[bucketCnt] = totalBucketCnt - curCnt;
        std::fill(fakeWeights.begin() + curCnt,
            fakeWeights.end(),
            (alpha / bucketCnt) / (1.0 / delta - 1.0));
        std::fill(realWeights.begin() + curCnt,
            realWeights.end(),
            0);
        for (uint32_t j = 0; j < replicaCnt[bucketCnt]; ++j)
            bucketIDList[curCnt + j] = makeBucketID(level, bucketCnt, j);
        std::partial_sum(replicaCnt.begin(), replicaCnt.end(), replicaCnt.begin());

        realDistributionLists[level] = std::discrete_distribution<uint32_t>(realWeights.begin(), realWeights.end());
        fakeDistributionLists[level] = std::discrete_distribution<uint32_t>(fakeWeights.begin(), fakeWeights.end());

        return replicaCnt;
    }
};
