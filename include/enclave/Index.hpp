#pragma once

#include "sgx_tcrypto.h"
#include <stdint.h>
#include <math.h> 
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
    int ocall_getBinData(const uint8_t* lvlCntPairs_, size_t size_,
        uint8_t** individualEncBins, uint32_t encBinSize,
        uint8_t** mergedBins, size_t* len,
        uint32_t binCnt,
        double storageOverhead);
    void ocall_insertBins(uint8_t lvl, uint32_t encBinSize,
        const uint8_t* individualEncBins,
        const uint8_t* i_j_, size_t size_);
    void ocall_delete(void* addr);
    void ocall_timing(uint8_t start, const char* msg);
    void ocall_setSeed(uint32_t seed);
};

class BaseIndex
{
public:
    virtual ~BaseIndex() {}
    virtual int insert(uint8_t* data) = 0;
    virtual int setup(uint8_t* data) = 0;
    virtual int pushDown() = 0;
    virtual size_t getBucketSize() const = 0;
    virtual uint32_t getBucketCapacity() const = 0;
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

    // differentially oblivious parameters
    const uint32_t lambda;
    const double epsilon;
    const double delta;

    // frequency smoothing
    const double storageOverhead;

    // parralelization 
    const uint32_t threadNum;

    // random seed 
    const uint32_t seed;

    const Datum<Key> dummyDatum;


    uint32_t binCapacity;
    uint32_t binSize;
    size_t encBinSize;

    uint8_t* secureBuf;

    std::vector<uint8_t> D;
    std::vector<std::vector<uint32_t>> replicaCntLists;
    std::vector<uint32_t> realCntList;
    std::mt19937 shuffleGen;

    Geom geom;

    uint32_t createReplicas(const std::vector<uint32_t>& replicaCnt,
        std::vector<std::pair<uint32_t, uint32_t>>& i_j,
        uint32_t originalBinCnt,
        uint8_t* individualEncBins)
    {
        static const double alpha = 1.0 / (storageOverhead - 1.0);
        static const double delta = (storageOverhead - 1.0) / storageOverhead;
        static const std::vector<uint8_t> encDummyBucket = [](const sgx_aes_ctr_128bit_key_t* table_key,
            uint32_t bucketCapacity, size_t encBucketSize, uint32_t binSize) {
                std::vector<Datum<Key>> dummyBin(bucketCapacity, Datum<Key>::getDummy());
                std::vector<uint8_t> ret(encBinSize);
                genNonce(ret.data(), SGX_AESGCM_IV_SIZE);
                sgx_rijndael128GCM_encrypt(table_key,
                    (uint8_t*)dummyBin.data(),
                    binSize,
                    ret.data() + SGX_AESGCM_IV_SIZE,
                    ret.data(), SGX_AESGCM_IV_SIZE,
                    NULL, 0,
                    (sgx_aes_gcm_128bit_tag_t*)(ret.data() + binSize + SGX_AESGCM_IV_SIZE)
                );
                return ret;
        }(this->table_key, this->bucketCapacity, this->getEncBucketSize, this->binSize);
        assert(encDummyBucket.size() == encBinSize);
        uint32_t totalBinCnt = std::accumulate(replicaCnt.begin(), replicaCnt.end(), 0);
        uint32_t binCnt = std::min(originalBinCnt, (uint32_t)replicaCnt.size());
        uint32_t curBinCnt = binCnt;
        uint8_t* targetEncBin = individualEncBins + (size_t)binCnt * encBinSize;
        uint8_t* sourceEncBin = individualEncBins;
        i_j.resize(totalBinCnt);
        for (uint32_t i = 0; i < binCnt; ++i)
        {
            i_j[i] = { i, 0 };
            for (uint32_t j = 1; j < replicaCnt[i]; ++j)
            {
                i_j[curBinCnt++] = { i, j };
                targetEncBin = std::copy(sourceEncBin, sourceEncBin + encBinSize, targetEncBin);
            }
            sourceEncBin += encBinSize;
        }

        assert(curBinCnt <= totalBinCnt);

        for (uint32_t i = curBinCnt; i < totalBinCnt; ++i)
        {
            i_j[i] = { binCnt, i - curBinCnt };
            targetEncBin = std::copy(encDummyBucket.data(), encDummyBucket.data() + encBinSize, targetEncBin);
        }
        return totalBinCnt;
    }

    void unBuild(const uint8_t& idx, std::vector<uint32_t>& replicas,
        uint32_t& pendingBinCnt, uint32_t& totalCnt, bool pushDown = false)
    {
        if (replicaCntLists[idx].empty())
        {
            if (pushDown)
                replicas.push_back(REPLICA_DEL);
            return;
        }
        replicas.insert(replicas.end(), replicaCntLists[idx].begin(), replicaCntLists[idx].end());
        replicas.push_back(REPLICA_DEL);
        pendingBinCnt += replicaCntLists[idx].size() - 1;
        totalCnt += realCntList[idx];
        replicaCntLists[idx].clear();
        realCntList[idx] = 0;
    }

    int allocateBins(uint8_t* mergedBins, uint8_t* individualEncBins,
        uint8_t& idx, uint32_t& binCnt, const uint32_t& totalCnt)
    {
        /**
         *  allocate data into bins
         */
        std::vector<uint32_t> binLoads(binCnt);

        // generate random bin loads 
        for (auto& l : binLoads)
            l = geom.sample();

        // compute the differentially private prefix sums
        uint32_t additiveError;
        std::vector<uint32_t> prefixSum = DPPrefixSum(binLoads, epsilon, seed,
            binCnt, binCapacity, additiveError);
        additiveError = upperBound<Datum<Key>>(additiveError);
        assert(binCapacity >= additiveError);

        std::vector<Datum<Key>> curBin(binCapacity);
        uint32_t binTagSize = binCnt * sizeof(BinTag<Key>);
        std::vector<uint8_t> binTags(binTagSize);
        BinTag<Key>* binTags_ = (BinTag<Key>*)binTags.data();

        uint32_t accCnt = additiveError;
        uint32_t curCnt = std::min(accCnt, totalCnt);
        uint32_t accLoad = 0;
        uint8_t* curMergedBin = mergedBins;
        uint8_t* curEncBin = individualEncBins;
        uint8_t ctr[16] = { 0 };
        Datum<Key>* secureDatum = (Datum<Key>*)secureBuf;
        std::fill(secureDatum, secureDatum + 3 * binCapacity + 2 * additiveError, dummyDatum);
        // put the first *s* data into the buffer 
        sgx_aes_ctr_decrypt(table_key, curMergedBin, datumSize * curCnt,
            ctr, 128, secureBuf);
        curMergedBin += datumSize * curCnt;
        uint32_t bufferLoad = curCnt;
        for (uint32_t i = 0; i < binCnt; ++i)
        {
            const uint32_t curLoad = accLoad + binLoads[i] > totalCnt ? totalCnt - accLoad : binLoads[i];
            curCnt = accCnt + prefixSum[i + 1] - prefixSum[i] > totalCnt ? totalCnt - accCnt : prefixSum[i + 1] - prefixSum[i];
            bufferLoad += curCnt;
            bufferLoad -= curLoad;
            if (curCnt > 0)
            {
                sgx_aes_ctr_decrypt(table_key,
                    curMergedBin, datumSize * curCnt,
                    ctr, 128,
                    (uint8_t*)(secureDatum + (binCapacity + 2 * additiveError))
                );
                accCnt += curCnt;
                curMergedBin += datumSize * curCnt;
            }
            accLoad += curLoad;

            // std::sort(secureDatum,
            //     secureDatum + binCapacity + 2 * additiveError + curCnt); // for comparing performance
            bitonicSort<Datum<Key>>(secureDatum, 0, binCapacity + 2 * additiveError + curCnt, true);

            assert(secureDatum[0].getDat() != 0);

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
            // for (uint32_t j = 0;j < 3; ++j)
            //     t_print("%d ", curBin[j].getDat());
            // binTags_[i] = DPInteriorPoint<Datum<Key>>(curBin, curLoad, epsilon, seed).getDat();
            if (curLoad > 0)
                binTags_[i] = { curBin[0].getDat(), curBin[curLoad - 1].getDat() };
            else
                binTags_[i] = { curBin[0].getDat(), curBin[0].getDat() };

            genNonce(curEncBin, SGX_AESGCM_IV_SIZE);
            sgx_rijndael128GCM_encrypt(table_key,
                (uint8_t*)curBin.data(), binSize,
                curEncBin + SGX_AESGCM_IV_SIZE,
                curEncBin, SGX_AESGCM_IV_SIZE,
                NULL, 0,
                (sgx_aes_gcm_128bit_tag_t*)(curEncBin + SGX_AESGCM_IV_SIZE + binSize));

            curEncBin += encBinSize;
        }

        realCntList[idx] = totalCnt;
        uint64_t sizeLevel = ((uint64_t)(binTagSize) << 32ll) | (uint64_t)idx;

        // write binTags to the client
        if (write_to_session_peer(
            sslSession, &sizeLevel, sizeof(sizeLevel)) != 0)
        {
            t_print(SERVER_PROXY " Write to client failed\n");
            return -1;
        }
        if (write_to_session_peer(
            sslSession, binTags.data(), binTagSize) != 0)
        {
            t_print(SERVER_PROXY " Write to client failed\n");
            return -1;
        }

        // read replicaCnt from the client
        size_t bytesToRead;
        if (read_from_session_peer(
            sslSession, &bytesToRead, sizeof(bytesToRead)) != sizeof(bytesToRead))
        {
            t_print(SERVER_PROXY " Read from client failed\n");
            return -1;
        }
        std::vector<uint32_t> replicaCnt(bytesToRead / sizeof(uint32_t));
        if (read_from_session_peer(
            sslSession, replicaCnt.data(), bytesToRead) != bytesToRead)
        {
            t_print(SERVER_PROXY " Read from client failed\n");
            return -1;
        }
        for (auto i = replicaCnt.rbegin(); i != std::prev(replicaCnt.rend()); ++i)
            *i -= *std::next(i);
        std::vector<std::pair<uint32_t, uint32_t>> i_j;
        uint32_t totalBinCnt = createReplicas(replicaCnt, i_j, binCnt, individualEncBins);
        replicaCntLists[idx] = std::move(replicaCnt);

        std::vector<uint32_t> shuffleKeys(totalBinCnt);
        std::iota(shuffleKeys.begin(), shuffleKeys.end(), 0);
        std::shuffle(shuffleKeys.begin(), shuffleKeys.end(), shuffleGen);

        bitonicSortWithBinsParallel(shuffleKeys.data(),
            totalBinCnt,
            table_key,
            i_j.data(),
            individualEncBins,
            encBinSize,
            threadNum
        );

        memcpy(mergedBins, (uint8_t*)i_j.data(), sizeof(i_j[0]) * i_j.size());
        ocall_insertBins(idx, encBinSize, individualEncBins,
            (uint8_t*)mergedBins, sizeof(i_j[0]) * i_j.size());
        return 0;
    }

    void retrievePendingBins(uint8_t& idx,
        uint32_t& pendingBinCnt,
        uint32_t& binCnt,
        uint32_t& totalCnt,
        uint8_t*& data,
        uint8_t*& individualEncBins,
        uint8_t*& mergedBins)
    {
        std::vector<uint32_t> replicas;
        unBuild(idx, replicas, pendingBinCnt, totalCnt);
        ++D[idx];

        while (D[idx] == D[idx + 1])
        {
            ++D[idx + 1];
            unBuild(idx, replicas, pendingBinCnt, totalCnt);
            D[idx] = idx;
            ++idx;
        }
        unBuild(idx, replicas, pendingBinCnt, totalCnt);

        size_t len = 0;

        binCnt = ceil(2.0 * totalCnt / (binCapacity * (1 - pow(log(lambda), -2))));

        ocall_getBinData((uint8_t*)replicas.data(),
            replicas.size() * sizeof(uint32_t),
            &individualEncBins, encBinSize,
            &mergedBins, &len,
            binCnt,
            storageOverhead);

        // for convenience, we append the newly arrived bin to the end of the buffer 
        genNonce(individualEncBins + len, SGX_AESGCM_IV_SIZE);
        sgx_rijndael128GCM_encrypt(table_key,
            data, binSize,
            individualEncBins + len + SGX_AESGCM_IV_SIZE,
            individualEncBins + len, SGX_AESGCM_IV_SIZE,
            NULL, 0,
            (sgx_aes_gcm_128bit_tag_t*)(individualEncBins + len + encBinSize - SGX_AESGCM_MAC_SIZE));
        len += encBinSize;
        ++pendingBinCnt;
    }

    void mergeBins(uint8_t* individualEncBins, uint8_t* mergedBins, const uint32_t pendingBinCnt)
    {
        // decrypt those bins and re-encrypt them as ctr mode 
        uint8_t ctr[16] = { 0 };
        uint8_t* curEncBin = individualEncBins;
        uint8_t* curMergedBin = mergedBins;

        assert(encBinSize - SGX_AESGCM_IV_SIZE - SGX_AESGCM_MAC_SIZE == binSize);
        for (uint32_t i = 0; i < pendingBinCnt; ++i)
        {
            sgx_rijndael128GCM_decrypt(table_key,
                curEncBin + SGX_AESGCM_IV_SIZE, binSize,
                secureBuf,
                curEncBin, SGX_AESGCM_IV_SIZE,
                NULL, 0,
                (sgx_aes_gcm_128bit_tag_t*)(curEncBin + encBinSize - SGX_AESGCM_MAC_SIZE));

            sgx_aes_ctr_encrypt(table_key,
                secureBuf, binSize,
                ctr, 128,
                curMergedBin);
            curEncBin += encBinSize;
            curMergedBin += binSize;
        }

        // obliviously sort those data 
        OSort<Datum<Key>>(table_key, (Datum<Key>*)secureBuf, DEFAULT_SECURE_BUFFER_SIZE,
            (Datum<Key>*)mergedBins, pendingBinCnt * binCapacity, threadNum);
    }

public:
    virtual size_t getBucketSize() const override
    {
        return binSize;
    }

    virtual uint32_t getBucketCapacity() const override
    {
        return binCapacity;
    }

    Index(const sgx_aes_ctr_128bit_key_t* table_key_,
        SSL*& sslSession,
        const uint32_t lambda, const double epsilon,
        const uint32_t binCapacity,
        const double delta, const uint32_t k,
        const double storageOverhead,
        const uint32_t threadNum, const uint32_t seed)
        : table_key(table_key_), sslSession(sslSession),
        lambda(lambda), epsilon(epsilon), delta(delta), storageOverhead(storageOverhead),
        threadNum(threadNum), seed(seed), dummyDatum(Datum<Key>::getDummy()), shuffleGen(seed)
    {
        // Compute bin capacity according to convolutional bound 
        // auto toEvenUB = [](const int x) { return x + (x & 1); };
        // uint32_t binCapacityUB = toEvenUB(std::floor(fastPower<double>(std::log(lambda), 5) / epsilon));
        // uint32_t binCapacityLB = 2;
        // while (binCapacityLB + 2 < binCapacityUB)
        // {
        //     binCapacity = toEvenUB((binCapacityLB + binCapacityUB) / 2);
        //     const std::vector<double> pmf = geomConv(binCapacity, epsilon, MIN_B);
        //     if (std::accumulate(pmf.begin(), pmf.begin() + binCapacity, 0.0) > delta)
        //         binCapacityLB = binCapacity + 2;
        //     else
        //         binCapacityUB = binCapacity;
        // }
        this->binCapacity = upperBound<Datum<Key>>(binCapacity);
        t_print("lambda = %d, epsilon = %f, delta = %f, k = %d, storageOverhead = %f, threadNum = %d, seed = %d\n",
            lambda, epsilon, delta, k, storageOverhead, threadNum, seed);
        // this->binCapacity = upperBound<Datum<Key>>(binCapacityUB);
        // this->binCapacity = std::max(binCapacity, 256u);
        binSize = this->binCapacity * datumSize;
        encBinSize = binSize + SGX_AESGCM_IV_SIZE + SGX_AESGCM_MAC_SIZE;

        replicaCntLists.resize(k + 1);
        realCntList.resize(k + 1);
        D.resize(k + 1);
        std::iota(D.begin(), D.end(), 0);
        D[k] = std::numeric_limits<uint8_t>::max();

        geom = Geom(binCapacity, epsilon, seed);
        genNonce(NULL, 0, seed); // initialize seed 

        secureBuf = new uint8_t[DEFAULT_SECURE_BUFFER_SIZE];
    }

    virtual int insert(uint8_t* data) override
    {
        uint8_t idx = 0;
        uint32_t pendingBinCnt = 0;
        uint32_t binCnt = 0;
        uint32_t totalCnt = binCapacity;
        uint8_t* individualEncBins = NULL;
        uint8_t* mergedBins = NULL;
        retrievePendingBins(idx,
            pendingBinCnt,
            binCnt, totalCnt,
            data,
            individualEncBins, mergedBins);


        mergeBins(individualEncBins, mergedBins, pendingBinCnt);

        int ret = allocateBins(mergedBins, individualEncBins, idx, binCnt, totalCnt);

        // ocall_delete(individualEncBins);
        // ocall_delete(mergedBins);
        return ret;
    }

    // setup from empty 
    virtual int setup(uint8_t* data) override
    {
        uint8_t idx = D.size() - 1;
        uint32_t binCnt = 0;
        uint32_t receivedCnt = 0;
        uint32_t totalCnt = *(uint32_t*)data;
        uint8_t* individualEncBins = NULL;
        uint8_t* mergedBins = NULL;

        size_t len = 0;
        binCnt = ceil(2.0 * totalCnt / (binCapacity * (1 - pow(log(lambda), -2))));

        std::vector<uint8_t> replicas;
        ocall_getBinData(replicas.data(),
            0,
            &individualEncBins, encBinSize,
            &mergedBins, &len,
            binCnt,
            storageOverhead);

        uint8_t ctr[16] = { 0 };
        uint8_t* curEncBin = mergedBins;
        uint32_t curCnt = std::min(binCapacity, totalCnt - receivedCnt);
        receivedCnt += curCnt;
        sgx_aes_ctr_encrypt(table_key,
            data + sizeof(uint32_t),
            sizeof(Datum<Key>) * curCnt,
            ctr, 128,
            curEncBin);
        curEncBin += sizeof(Datum<Key>) * curCnt;
        std::vector<uint8_t> buffer;

        while (receivedCnt < totalCnt)
        {
            curCnt = std::min(binCapacity, totalCnt - receivedCnt);
            receivedCnt += curCnt;
            size_t bytesToRead = 0;
            if (read_from_session_peer(
                sslSession, &bytesToRead, sizeof(bytesToRead)) != sizeof(bytesToRead))
            {
                t_print(SERVER_PROXY " Read from client failed\n");
                return -1;
            }
            buffer.resize(bytesToRead);
            if (read_from_session_peer(
                sslSession, buffer.data(), bytesToRead) != bytesToRead)
            {
                t_print(SERVER_PROXY " Read from client failed\n");
                return -1;
            }
            sgx_aes_ctr_encrypt(table_key,
                buffer.data() + sizeof(OP) + sizeof(uint32_t),
                sizeof(Datum<Key>) * curCnt,
                ctr, 128,
                curEncBin);
            curEncBin += sizeof(Datum<Key>) * curCnt;
        }
        OSort<Datum<Key>>(table_key,
            (Datum<Key>*)secureBuf, DEFAULT_SECURE_BUFFER_SIZE,
            (Datum<Key>*)mergedBins, totalCnt,
            threadNum);

        int ret = allocateBins(mergedBins, individualEncBins, idx, binCnt, totalCnt);

        write_to_session_peer(sslSession, &ret, sizeof(ret));
        // ocall_delete(individualEncBins);
        // ocall_delete(mergedBins);
        return ret;
    }

    // merge all bin lists into one bin list
    virtual int pushDown() override
    {
        uint8_t idx = 0;
        uint32_t pendingBinCnt = 0;
        uint32_t binCnt = 0;
        uint32_t totalCnt = 0;
        uint8_t* individualEncBins = NULL;
        uint8_t* mergedBins = NULL;

        std::vector<uint32_t> replicas;

        for (uint8_t i = 0; i < D.size(); ++i)
        {
            D[i] = i;
            unBuild(i, replicas, pendingBinCnt, totalCnt, true);
        }

        D.back() = std::numeric_limits<uint8_t>::max();
        idx = D.size() - 1;

        size_t len = 0;
        binCnt = ceil(2.0 * totalCnt / (binCapacity * (1 - pow(log(lambda), -2))));


        ocall_getBinData((uint8_t*)replicas.data(),
            replicas.size() * sizeof(uint32_t),
            &individualEncBins, encBinSize,
            &mergedBins, &len,
            binCnt,
            storageOverhead);

        mergeBins(individualEncBins, mergedBins, pendingBinCnt);

        int ret = allocateBins(mergedBins, individualEncBins, idx, binCnt, totalCnt);

        write_to_session_peer(sslSession, &ret, sizeof(ret));
        // ocall_delete(individualEncBins);
        // ocall_delete(mergedBins);
        return ret;
    }


    virtual ~Index()
    {
        delete[] secureBuf;
    }
};
