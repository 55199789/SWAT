#pragma once

#include "sgx_tcrypto.h"
#include <math.h> 
#include <stdint.h>
#include <algorithm>
#include <deque>
#include <limits>
#include <random>
#include <typeinfo>
#include <vector>
#include "common/common.h"
#include "common/defs.h"
#include "common/openssl_utility.h"
#include "enclave/Dist.hpp"
#include "enclave/DPInteriorPoint.hpp"
#include "enclave/DPPrefixSum.hpp"
#include "enclave/ObliviousSort.hpp"

extern "C"
{
    void ocall_resizeMemoryPool(uint8_t** ptr, uint32_t poolIndex, size_t size);
    void ocall_getBucketData(const uint32_t* replicaCnts, size_t rCnt,
        const uint32_t* bucketCnts, size_t bCnt, uint8_t** encBucket,
        uint32_t encBucketSize);
    void ocall_insertBuckets(uint8_t lvl, uint32_t encBucketSize,
        const uint8_t* individualEncBuckets,
        const uint8_t* i_j_, size_t size_);
    void ocall_delete(void* addr);
    void ocall_timing(uint8_t start, const char* msg);
    void ocall_setSeed(uint32_t seed);
};

class BaseIndex
{
public:
    virtual ~BaseIndex() {}
    virtual int insert(uint8_t* encData, uint32_t initBucketCnt) = 0;
    virtual size_t getBucketSize() const = 0;
    virtual size_t getEncBucketSize() const = 0;
    virtual size_t getBucketCapacity() const = 0;
    virtual size_t getDatumSize() const = 0;
    virtual void osort(uint8_t* encData, uint32_t cnt) = 0;
};

/**
 * @brief Differentially oblivious index structure
 * of Key type
 *
 * @tparam Key, integer or float numbers
 */
template <class Key>
class Index : public BaseIndex
{
private:
    static constexpr uint32_t datumSize = sizeof(Datum<Key>);

    const sgx_aes_ctr_128bit_key_t* table_key;
    SSL* sslSession;
    std::vector<uint8_t> secureBuf;

    // differentially oblivious parameters
    const uint32_t lambda;
    const double epsilon;
    const double delta;

    // bucketization 
    const uint32_t bucketCapacity;
    const uint32_t bucketSize;
    const uint32_t encBucketSize;
    std::mt19937 shuffleGen;

    // frequency smoothing
    const double storageOverhead;

    // parralelization 
    const uint32_t threadNum;

    // random seed 
    const uint32_t seed;
    const Datum<Key> dummyDatum;

    // dynamization
    std::vector<uint8_t> D;
    uint32_t binCapacity;
    uint32_t binSize;
    uint32_t encBinSize;
    Geom geom;

    // frequency smoothing
    std::vector<std::vector<uint32_t>> replicaCntLists;

    // auxiliary data 
    double factor;
    std::mutex samplingMutex;

    void unBuild(const uint8_t& idx, std::vector<uint32_t>& replicas,
        std::vector<uint32_t>& pendingBucketCnts,
        std::vector<uint32_t>& binCnts)
    {
        if (replicaCntLists[idx].empty())
            return;
        replicas.insert(replicas.end(), replicaCntLists[idx].begin(), replicaCntLists[idx].end());
        pendingBucketCnts.push_back(replicaCntLists[idx].size() - 1);
        replicaCntLists[idx].clear();
        binCnts.push_back(
            (uint32_t)std::ceil(
                factor * pendingBucketCnts.back()
            )
        );
    }


    uint32_t createReplicas(const std::vector<uint32_t>& replicaCnts,
        std::vector<std::pair<uint32_t, uint32_t>>& i_j,
        uint8_t* encBuckets)
    {
        static const double alpha = 1.0 / (storageOverhead - 1.0);
        static const double delta = (storageOverhead - 1.0) / storageOverhead;
        static const std::vector<uint8_t> encDummyBucket = [](const sgx_aes_ctr_128bit_key_t* table_key,
            uint32_t bucketCapacity, size_t encBucketSize, uint32_t bucketSize) {
                std::vector<Datum<Key>> dummyBucket(bucketCapacity, Datum<Key>::getDummy());
                std::vector<uint8_t> ret(encBucketSize);
                genNonce(ret.data(), SGX_AESGCM_IV_SIZE);
                sgx_rijndael128GCM_encrypt(table_key,
                    (uint8_t*)dummyBucket.data(),
                    bucketSize,
                    ret.data() + SGX_AESGCM_IV_SIZE,
                    ret.data(), SGX_AESGCM_IV_SIZE,
                    NULL, 0,
                    (sgx_aes_gcm_128bit_tag_t*)(ret.data() + bucketSize + SGX_AESGCM_IV_SIZE)
                );
                return ret;
        }(this->table_key, this->bucketCapacity, encBucketSize, this->bucketSize);

        if (replicaCnts.empty())
            return 0;

        const uint32_t duplicatedBucketCnt = std::accumulate(replicaCnts.begin(), replicaCnts.end(), 0);
        const uint32_t totalBucketCnt = replicaCnts.size() - 1;
        uint32_t curBucketCnt = totalBucketCnt;
        uint8_t* targetEncBucket = encBuckets + (size_t)totalBucketCnt * encBucketSize;
        uint8_t* sourceEncBucket = encBuckets;
        i_j.resize(duplicatedBucketCnt);
        for (uint32_t i = 0; i < totalBucketCnt; ++i)
        {
            i_j[i] = { i, 0 };
            for (uint32_t j = 1; j < replicaCnts[i]; ++j)
            {
                i_j[curBucketCnt++] = { i, j };
                targetEncBucket = std::copy(sourceEncBucket, sourceEncBucket + encBucketSize, targetEncBucket);
            }
            sourceEncBucket += encBucketSize;
        }

        for (uint32_t i = curBucketCnt; i < duplicatedBucketCnt; ++i)
        {
            i_j[i] = { totalBucketCnt, i - curBucketCnt };
            targetEncBucket = std::copy(encDummyBucket.data(), encDummyBucket.data() + encBucketSize, targetEncBucket);
        }
        return duplicatedBucketCnt;
    }

    uint32_t DOAllocate(uint8_t* encData, uint8_t* encBinBuffer, uint32_t dataCnt,
        std::vector<uint32_t>& prefixSums,
        std::vector<std::pair<Key, int>>& interiorPoints,
        uint32_t& additiveError
    )
    {
        const uint32_t binCnt = (uint32_t)std::ceil(
            factor * (dataCnt / bucketCapacity)
        );
        std::vector<uint32_t> binLoads(binCnt);
        std::vector<Datum<Key>> curBin(binCapacity);
        sgx_status_t ret;

        // generate random bin loads 
        {
#ifdef _OPENMP // for thread safety
            const std::lock_guard<std::mutex> lock(samplingMutex);
#endif 
            uint32_t totalLoad = 0;
            for (uint32_t i = 0;i < binCnt; ++i)
            {
                binLoads[i] = std::min(geom.sample(), (uint32_t)(dataCnt - totalLoad));
                totalLoad += binLoads[i];
            }
            prefixSums = DPPrefixSum(binLoads, epsilon, seed,
                binCapacity, additiveError);
        }
        interiorPoints.resize(binCnt);

        // compute the differentially private prefix sums
        additiveError = upperBound(datumSize, additiveError);

        uint32_t accCnt = additiveError;
        uint32_t curCnt = std::min(accCnt, dataCnt);
        uint32_t accLoad = 0;
        uint8_t* curEncData = encData;
        uint8_t* curEncBin = encBinBuffer;
        uint8_t ctr[16] = { 0 };
        std::vector<Datum<Key>> secureDatumVec(3 * binCapacity + 2 * additiveError, dummyDatum);
        Datum<Key>* secureDatum = secureDatumVec.data();
        // put the first *s* data into the buffer 
        ret = sgx_aes_ctr_decrypt(table_key, curEncData, datumSize * curCnt,
            ctr, 128, (uint8_t*)secureDatum);
        if (ret != SGX_SUCCESS)
        {
            t_print("DOAllocate: sgx_aes_ctr_decrypt failed: %d\n", ret);
            return ret;
        }
        curEncData += datumSize * curCnt;
        for (uint32_t i = 0; i < binCnt; ++i)
        {
            const uint32_t curLoad = accLoad + binLoads[i] > dataCnt ? dataCnt - accLoad : binLoads[i];
            curCnt = accCnt + prefixSums[i + 1] - prefixSums[i] > dataCnt ? dataCnt - accCnt : prefixSums[i + 1] - prefixSums[i];
            if (curCnt > 0)
            {
                ret = sgx_aes_ctr_decrypt(table_key,
                    curEncData, datumSize * curCnt,
                    ctr, 128,
                    (uint8_t*)(secureDatum + (binCapacity + 2 * additiveError))
                );
                if (ret != SGX_SUCCESS)
                {
                    t_print("DOAllocate: sgx_aes_ctr_decrypt %d failed: %d\n", i, ret);
                    return ret;
                }
                accCnt += curCnt;
                curEncData += datumSize * curCnt;
            }
            accLoad += curLoad;

            // std::sort(secureDatum,
            //     secureDatum + binCapacity + 2 * additiveError + curCnt); // for comparing performance
            bitonicSort<Datum<Key>>(secureDatum, 0, binCapacity + 2 * additiveError + curCnt, true);

            for (uint32_t j = 0; j < binCapacity; ++j)
                if (j < curLoad)
                {
                    curBin[j] = secureDatum[j];
                    secureDatum[j] = dummyDatum;
                }
                else
                {
                    curBin[j] = dummyDatum;
                    secureDatum[j] = secureDatum[j];
                }
            interiorPoints[i] = {
                DPInteriorPoint<Datum<Key>>(curBin, curLoad, epsilon, seed).getDat(),
                (int)i + 1
            };

            genNonce(curEncBin, SGX_AESGCM_IV_SIZE);
            ret = sgx_rijndael128GCM_encrypt(table_key,
                (uint8_t*)curBin.data(), binSize,
                curEncBin + SGX_AESGCM_IV_SIZE,
                curEncBin, SGX_AESGCM_IV_SIZE,
                NULL, 0,
                (sgx_aes_gcm_128bit_tag_t*)(curEncBin + SGX_AESGCM_IV_SIZE + binSize));
            if (ret != SGX_SUCCESS)
            {
                t_print("DOAllocate: sgx_rijndael128GCM_encrypt failed: %d\n", ret);
                return ret;
            }
            curEncBin += encBinSize;
        }
        return binCnt;
    }

    void DOMerge(uint8_t* leftArr, uint8_t* rightArr,
        uint8_t* leftEncBins, uint8_t* rightEncBins,
        uint32_t leftCnt, uint32_t rightCnt)
    {
        // allocate sorted arrays into bins 
        uint32_t leftAdditiveError, rightAdditiveError, additiveError;
        uint32_t leftBinCnt, rightBinCnt;
        std::vector<uint32_t> leftPrefixSums;
        std::vector<std::pair<Key, int>> leftInteriorPoints;
        std::vector<uint32_t> rightPrefixSums;
        std::vector<std::pair<Key, int>> rightInteriorPoints;
#ifdef _OPENMP 
#pragma omp taskgroup
        {
#pragma omp task default(shared)
            {
#endif 
                leftBinCnt = DOAllocate(leftArr, leftEncBins, leftCnt,
                    leftPrefixSums, leftInteriorPoints,
                    leftAdditiveError);
#ifdef _OPENMP
            }
#pragma omp task default(shared)
            {
#endif
                rightBinCnt = DOAllocate(rightArr, rightEncBins, rightCnt,
                    rightPrefixSums, rightInteriorPoints,
                    rightAdditiveError);
#ifdef _OPENMP 
            }
#pragma omp taskyield
        }
#endif 
        additiveError = std::max(leftAdditiveError, rightAdditiveError);
        // merge bins 
        sgx_status_t ret;
        std::vector<Datum<Key>> secureDatumVec(7 * (binCapacity + additiveError), dummyDatum);
        std::vector<std::pair<Key, int>> interiorPoints = leftInteriorPoints;
        Datum<Key>* secureDatum = secureDatumVec.data();
        int j0 = -1, j1 = -1;
        int k0 = 0, k1 = 0;
        uint32_t binCnt = leftBinCnt + rightBinCnt;
        uint32_t cnt = 0;
        uint8_t* leftEncBin = leftEncBins;
        uint8_t* rightEncBin = rightEncBins;
        uint8_t* mergedArr = leftArr;
        uint8_t ctr[16] = { 0 };

        for (auto& i : rightInteriorPoints)
            i.second = -i.second;
        interiorPoints.insert(interiorPoints.end(),
            rightInteriorPoints.begin(), rightInteriorPoints.end());

        std::sort(interiorPoints.begin(), interiorPoints.end());

        ret = sgx_rijndael128GCM_decrypt(
            table_key,
            leftEncBin + SGX_AESGCM_IV_SIZE, binSize,
            (uint8_t*)secureDatum,
            leftEncBin, SGX_AESGCM_IV_SIZE,
            NULL, 0,
            (sgx_aes_gcm_128bit_tag_t*)(leftEncBin + SGX_AESGCM_IV_SIZE + binSize)
        );
        if (ret != SGX_SUCCESS)
        {
            t_print("DOMerge: sgx_rijndael128GCM_decrypt left failed: %d\n", ret);
            return;
        }
        leftEncBin += encBinSize;

        ret = sgx_rijndael128GCM_decrypt(
            table_key,
            rightEncBin + SGX_AESGCM_IV_SIZE, binSize,
            (uint8_t*)(secureDatum + binCapacity),
            rightEncBin, SGX_AESGCM_IV_SIZE,
            NULL, 0,
            (sgx_aes_gcm_128bit_tag_t*)(rightEncBin + SGX_AESGCM_IV_SIZE + binSize)
        );
        if (ret != SGX_SUCCESS)
        {
            t_print("DOMerge: sgx_rijndael128GCM_decrypt right failed: %d\n", ret);
            return;
        }
        rightEncBin += encBinSize;

        for (uint32_t i = 0; i < binCnt; ++i)
        {
            bool left = interiorPoints[i].second > 0;
            if (left && j0 + 2 < (int)leftBinCnt) // left array
            {
                ++j0;
                ret = sgx_rijndael128GCM_decrypt(
                    table_key,
                    leftEncBin + SGX_AESGCM_IV_SIZE, binSize,
                    (uint8_t*)(secureDatum + 6 * (binCapacity + additiveError)),
                    leftEncBin, SGX_AESGCM_IV_SIZE,
                    NULL, 0,
                    (sgx_aes_gcm_128bit_tag_t*)(leftEncBin + SGX_AESGCM_IV_SIZE + binSize)
                );
                if (ret != SGX_SUCCESS)
                {
                    t_print("DOMerge: sgx_rijndael128GCM_decrypt left %d failed: %d\n", j0, ret);
                    return;
                }

                leftEncBin += encBinSize;

            }
            else if (j1 + 2 < (int)rightBinCnt) // right array
            {
                ++j1;
                ret = sgx_rijndael128GCM_decrypt(
                    table_key,
                    rightEncBin + SGX_AESGCM_IV_SIZE, binSize,
                    (uint8_t*)(secureDatum + 6 * (binCapacity + additiveError)),
                    rightEncBin, SGX_AESGCM_IV_SIZE,
                    NULL, 0,
                    (sgx_aes_gcm_128bit_tag_t*)(rightEncBin + SGX_AESGCM_IV_SIZE + binSize)
                );
                if (ret != SGX_SUCCESS)
                {
                    t_print("DOMerge: sgx_rijndael128GCM_decrypt right %d failed: %d\n", j1, ret);
                    return;
                }

                rightEncBin += encBinSize;
            }

            bitonicSort<Datum<Key>>(secureDatum, 0, 7 * binCapacity + 6 * additiveError, true);

            // determine safe bins 
            k0 = (int)leftBinCnt - 2 > j0 + 1 ? j0 + 1 : (int)leftBinCnt - 2;
            while (k0 > -1 && rightInteriorPoints[j1 + 1].first < leftInteriorPoints[k0 + 1].first)
                --k0;
            k1 = (int)rightBinCnt - 2 > j1 + 1 ? j1 + 1 : (int)rightBinCnt - 2;
            while (k1 > -1 && leftInteriorPoints[j0 + 1].first < rightInteriorPoints[k1 + 1].first)
                --k1;

            // put safe records into merged array
            uint32_t newCnt = leftPrefixSums[k0 + 1] + rightPrefixSums[k1 + 1] > 2 * additiveError ?
                leftPrefixSums[k0 + 1] + rightPrefixSums[k1 + 1] - 2 * additiveError :
                0;
            if (newCnt > cnt)
            {
                ret = sgx_aes_ctr_encrypt(
                    table_key,
                    (uint8_t*)secureDatum, (newCnt - cnt) * datumSize,
                    ctr, 128,
                    mergedArr
                );
                if (ret != SGX_SUCCESS)
                {
                    t_print("DOMerge: sgx_aes_ctr_encrypt %d failed: %d\n", i, ret);
                    return;
                }
                std::fill(secureDatum, secureDatum + (newCnt - cnt), dummyDatum);
                for (uint32_t j = 0; j < std::min(binCapacity, newCnt - cnt); ++j)
                {
                    std::swap(secureDatum[j], secureDatum[6 * (binCapacity + additiveError) + j]);
                }

                mergedArr += (newCnt - cnt) * datumSize;
                cnt = newCnt;
            }
        }
        if (cnt < leftCnt + rightCnt)
        {
            bitonicSort<Datum<Key>>(secureDatum, 0, 6 * (binCapacity + additiveError), true);
            ret = sgx_aes_ctr_encrypt(
                table_key,
                (uint8_t*)secureDatum, (leftCnt + rightCnt - cnt) * datumSize,
                ctr, 128,
                mergedArr
            );
            if (ret != SGX_SUCCESS)
            {
                t_print("DOMerge: sgx_aes_ctr_encrypt failed: %d\n", ret);
                return;
            }
        }
    }

    void KWayDOMerge(const std::vector<uint8_t*>& sortedArrays,
        const std::vector<uint8_t*>& binPtrs,
        std::vector<uint32_t>& pendingBucketCnts,
        uint32_t left, uint32_t right)
    {
        if (left == right)
            return;
        if (left + 1 == right)
        {
            DOMerge(sortedArrays[left], sortedArrays[right],
                binPtrs[left], binPtrs[right],
                pendingBucketCnts[left] * bucketCapacity,
                pendingBucketCnts[right] * bucketCapacity);
            pendingBucketCnts[left] += pendingBucketCnts[right];
            return;
        }

#ifdef _OPENMP 
#pragma omp taskgroup
        {
#pragma omp task default(shared)
#endif 
            KWayDOMerge(sortedArrays,
                binPtrs,
                pendingBucketCnts,
                left, (left + right) / 2);
#ifdef _OPENMP
#pragma omp task default(shared)
#endif
            KWayDOMerge(sortedArrays,
                binPtrs,
                pendingBucketCnts,
                (left + right) / 2 + 1, right);
#ifdef _OPENMP 
#pragma omp taskyield
        }
#endif 

        DOMerge(sortedArrays[left],
            sortedArrays[(left + right) / 2 + 1],
            binPtrs[left],
            binPtrs[(left + right) / 2 + 1],
            pendingBucketCnts[left] * bucketCapacity,
            pendingBucketCnts[(left + right) / 2 + 1] * bucketCapacity);
        pendingBucketCnts[left] += pendingBucketCnts[(left + right) / 2 + 1];
    }

public:

    Index(const sgx_aes_ctr_128bit_key_t* table_key_,
        SSL*& sslSession,
        const uint32_t lambda, const double epsilon, const double delta,
        const uint32_t bucketCapacity,
        const double storageOverhead,
        const uint32_t threadNum, const uint32_t seed,
        const uint32_t k)
        : table_key(table_key_), sslSession(sslSession),
        secureBuf(DEFAULT_SECURE_BUFFER_SIZE),
        lambda(lambda), epsilon(epsilon), delta(delta),
        bucketCapacity(bucketCapacity),
        bucketSize(bucketCapacity* datumSize),
        encBucketSize(bucketSize + SGX_AESGCM_IV_SIZE + SGX_AESGCM_MAC_SIZE),
        shuffleGen(seed),
        storageOverhead(storageOverhead),
        threadNum(threadNum), seed(seed),
        dummyDatum(Datum<Key>::getDummy()),
        D(k + 1),
        replicaCntLists(k)
    {
        t_print("Datum size = %u\n", datumSize);
        assert(getBucketSize() == upperBound(datumSize, getBucketSize()));
        // Compute bin capacity according to convolutional bound 
        auto toEvenUB = [](const int x) { return x + (x & 1); };
        uint32_t binCapacityUB = toEvenUB(std::floor(fastPower<double>(std::log(lambda), 5) / epsilon));
        uint32_t binCapacityLB = 2;
        while (binCapacityLB + 2 < binCapacityUB)
        {
            binCapacity = toEvenUB((binCapacityLB + binCapacityUB) / 2);
            uint32_t minB = std::ceil(2 * bucketCapacity / (binCapacity * (1 - pow(log(lambda), -2))));
            const std::vector<long double> pmf = Geom::geomConv(binCapacity, epsilon, minB);
            if (std::accumulate(pmf.begin(), pmf.begin() + bucketCapacity, 0.0) > delta)
                binCapacityLB = binCapacity + 2;
            else
                binCapacityUB = binCapacity;
        }
        binCapacity = upperBound(datumSize, binCapacity);
        binSize = binCapacity * datumSize;
        encBinSize = binSize + SGX_AESGCM_IV_SIZE + SGX_AESGCM_MAC_SIZE;

        std::iota(D.begin(), D.end(), 0);
        D[k] = std::numeric_limits<uint8_t>::max();
        geom = Geom(binCapacity, epsilon, seed);

        factor = 2.0 * bucketCapacity \
            / (binCapacity * (1.0 - pow(log(lambda), -2)));

        if (write_to_session_peer(
            sslSession, &binCapacity, sizeof(binCapacity)) != 0)
        {
            t_print(SERVER_PROXY " Write to client failed\n");
            assert(false);
        }
    }


    virtual void osort(uint8_t* encData, uint32_t cnt) override
    {
        OSort<Datum<Key>>(table_key,
            (Datum<Key>*)secureBuf.data(), DEFAULT_SECURE_BUFFER_SIZE,
            (Datum<Key>*)encData, cnt, threadNum);
    }

    virtual int insert(uint8_t* encData, uint32_t initBucketCnt) override
    {
        uint8_t idx = 0;
        std::vector<uint32_t> replicas;
        std::vector<uint32_t> pendingBucketCnts = { initBucketCnt };
        std::vector<uint32_t> binCnts = { (uint32_t)std::ceil(factor * initBucketCnt) };
        // unBuild(idx, replicas, pendingBucketCnts);
        ++D[idx];

        while (D[idx] == D[idx + 1])
        {
            ++D[idx + 1];
            unBuild(idx, replicas, pendingBucketCnts, binCnts);
            D[idx] = idx;
            ++idx;
        }
        unBuild(idx, replicas, pendingBucketCnts, binCnts);

        uint32_t totalBucketCnt =
            std::accumulate(pendingBucketCnts.begin(), pendingBucketCnts.end(), 0);
        ocall_getBucketData(
            replicas.data(), replicas.size(),
            pendingBucketCnts.data(), pendingBucketCnts.size(),
            &encData, encBucketSize
        );

        // Transform GCM ciphertexts to CTR ciphertexts
        uint8_t ctr[16] = { 0 };
#ifdef OBLIVIOUSMERGE
        ctr128_inc_for_bytes(ctr, bucketSize * initBucketCnt);
#endif 
        uint8_t* curGCMEncBucket = encData + bucketSize * initBucketCnt;
        uint8_t* curCTREncBucket = curGCMEncBucket;
        for (uint32_t i = 1;i < pendingBucketCnts.size();++i)
        {
#ifndef OBLIVIOUSMERGE
            memset(ctr, 0, sizeof(ctr));
#endif 
            for (uint32_t j = 0;j < pendingBucketCnts[i];++j)
            {
                // decrypt the bucket into the secure buffer
                sgx_status_t ret = sgx_rijndael128GCM_decrypt(
                    table_key,
                    curGCMEncBucket + SGX_AESGCM_IV_SIZE,
                    bucketSize,
                    (uint8_t*)secureBuf.data(),
                    curGCMEncBucket, SGX_AESGCM_IV_SIZE,
                    NULL, 0,
                    (const sgx_aes_gcm_128bit_tag_t*)(curGCMEncBucket + SGX_AESGCM_IV_SIZE + bucketSize)
                );
                if (ret != SGX_SUCCESS)
                {
                    t_print("sgx_rijndael128GCM_decrypt failed: %d\n", ret);
                    return ret;
                }

                // encrypt the bucket from the secure buffer to the untrusted memory
                ret = sgx_aes_ctr_encrypt(
                    table_key,
                    (uint8_t*)secureBuf.data(),
                    bucketSize,
                    ctr,
                    128,
                    curCTREncBucket
                );
                if (ret != SGX_SUCCESS)
                {
                    t_print("sgx_aes_ctr_encrypt failed: %d\n", ret);
                    return ret;
                }

                curCTREncBucket += bucketSize;
                curGCMEncBucket += encBucketSize;
            }
        }

        // merge those data 
        uint8_t* encBinData = NULL;
#ifndef OBLIVIOUSMERGE
        // Differentially oblivious merge sort
        ocall_resizeMemoryPool(&encBinData, 0,
            (size_t)std::accumulate(
                binCnts.begin(), binCnts.end(), 0
            ) * (size_t)encBinSize
        );
        std::vector<uint8_t*> sortedArrays;
        std::vector<uint8_t*> binPtrs;
        uint32_t accBucketCnts = 0;
        uint32_t accBinCnts = 0;
        for (uint32_t i = 0; i < binCnts.size(); ++i)
        {
            sortedArrays.push_back(encData + accBucketCnts * bucketSize);
            binPtrs.push_back(encBinData + accBinCnts * encBinSize);

            accBucketCnts += pendingBucketCnts[i];
            accBinCnts += binCnts[i];
        }

#ifdef _OPENMP 
#pragma omp parallel num_threads(threadNum) if (binCnts.size() > 2)
#pragma omp single
        {
#endif
            KWayDOMerge(sortedArrays, binPtrs,
                pendingBucketCnts,
                0, binCnts.size() - 1);
#ifdef _OPENMP
        }
#endif 
#endif
#ifdef OBLIVIOUSMERGE
        // Oblivious merge sort
        OSort<Datum<Key>>(table_key,
            (Datum<Key>*)secureBuf.data(), DEFAULT_SECURE_BUFFER_SIZE,
            (Datum<Key>*)encData, totalBucketCnt * bucketCapacity,
            threadNum);
#endif 
        ocall_resizeMemoryPool(&encBinData, 0,
            (size_t)std::ceil((size_t)totalBucketCnt * encBucketSize * storageOverhead)
        );

        // Bucketization
        uint32_t bucketTagSize = totalBucketCnt * sizeof(BucketTag<Key>);
        std::vector<uint8_t> bucketTags(bucketTagSize);
        BucketTag<Key>* bucketTags_ = (BucketTag<Key>*)bucketTags.data();
        curGCMEncBucket = encBinData;
        curCTREncBucket = encData;
        uint64_t sizeLevel = ((uint64_t)(bucketTagSize) << 32ll) | (uint64_t)idx;
        memset(ctr, 0, sizeof(ctr));
        Datum<Key>* curBucketData = (Datum<Key>*)secureBuf.data();
        for (uint32_t i = 0; i < totalBucketCnt; ++i)
        {
            // decrypt a bucket of records into the secure buffer 
            sgx_status_t ret = sgx_aes_ctr_decrypt(
                table_key,
                curCTREncBucket, bucketSize,
                ctr, 128,
                (uint8_t*)curBucketData
            );
            if (ret != SGX_SUCCESS)
            {
                t_print("insert: sgx_aes_ctr_decrypt %d failed: %d\n", i, ret);
                return ret;
            }

            bucketTags_[i] = {
                curBucketData[0].getDat(),
                curBucketData[bucketCapacity - 1].getDat()
            };

            // encrypt it into the untrusted memory
            genNonce(curGCMEncBucket, SGX_AESGCM_IV_SIZE);
            ret = sgx_rijndael128GCM_encrypt(
                table_key,
                (uint8_t*)curBucketData,
                bucketSize,
                (uint8_t*)curGCMEncBucket + SGX_AESGCM_IV_SIZE,
                curGCMEncBucket, SGX_AESGCM_IV_SIZE,
                NULL, 0,
                (sgx_aes_gcm_128bit_tag_t*)(curGCMEncBucket + SGX_AESGCM_IV_SIZE + bucketSize)
            );
            if (ret != SGX_SUCCESS)
            {
                t_print("insert: sgx_rijndael128GCM_encrypt %d failed: %d\n", i, ret);
                return ret;
            }

            curCTREncBucket += bucketSize;
            curGCMEncBucket += encBucketSize;
        }

        // write bucketTags to the client
        if (write_to_session_peer(
            sslSession, &sizeLevel, sizeof(sizeLevel)) != 0)
        {
            t_print(SERVER_PROXY " Write to client failed\n");
            return -1;
        }
        if (write_to_session_peer(
            sslSession, bucketTags.data(), bucketTagSize) != 0)
        {
            t_print(SERVER_PROXY " Write to client failed\n");
            return -1;
        }

        // read replicaCnts from the client
        size_t bytesToRead;
        if (read_from_session_peer(
            sslSession, &bytesToRead, sizeof(bytesToRead)) != sizeof(bytesToRead))
        {
            t_print(SERVER_PROXY " Read from client failed\n");
            return -1;
        }
        std::vector<uint32_t> replicaCnts(bytesToRead / sizeof(uint32_t));
        if (read_from_session_peer(
            sslSession, replicaCnts.data(), bytesToRead) != bytesToRead)
        {
            t_print(SERVER_PROXY " Read from client failed\n");
            return -1;
        }
        for (auto i = replicaCnts.rbegin(); i != std::prev(replicaCnts.rend()); ++i)
            *i -= *std::next(i);
        std::vector<std::pair<uint32_t, uint32_t>> i_j;
        const uint32_t duplicatedBucketCnt = createReplicas(replicaCnts, i_j, encBinData);
        replicaCntLists[idx] = std::move(replicaCnts);

        std::vector<uint32_t> shuffleKeys(duplicatedBucketCnt);
        std::iota(shuffleKeys.begin(), shuffleKeys.end(), 0);
        std::shuffle(shuffleKeys.begin(), shuffleKeys.end(), shuffleGen);
        bitonicSortWithBucketsParallel(
            secureBuf.data(),
            shuffleKeys.data(), duplicatedBucketCnt,
            table_key,
            i_j.data(),
            encBinData,
            encBucketSize,
            threadNum
        );


        memcpy(encData, (uint8_t*)i_j.data(), sizeof(i_j[0]) * i_j.size());
        ocall_insertBuckets(idx, encBucketSize, encBinData,
            (uint8_t*)encData, sizeof(i_j[0]) * i_j.size());

        int ret = 0;
        if (write_to_session_peer(
            sslSession, &ret, sizeof(ret)) != 0)
        {
            t_print(SERVER_PROXY " Write to client failed\n");
            return -1;
        }
        return 0;
    }

    virtual size_t getBucketSize() const override
    {
        return bucketSize;
    }

    virtual size_t getEncBucketSize() const override
    {
        return encBucketSize;
    }

    virtual size_t getBucketCapacity() const override
    {
        return bucketCapacity;
    }

    virtual size_t getDatumSize() const override
    {
        return datumSize;
    }

    virtual ~Index()
    {
    }
};
