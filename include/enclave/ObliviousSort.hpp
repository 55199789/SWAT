#pragma once 

#include <math.h>
#include <atomic>
#include <algorithm>
#include <thread>

void genNonce(uint8_t* nonce, const uint32_t len, const uint32_t seed = 0)
{
    // sgx_read_rand(nonce, len);
    static std::mt19937 gen(seed);
    static std::uniform_int_distribution<int> dist(0, 0x7fffffff);
    for (uint32_t i = 0; i < len / sizeof(int); ++i)
        ((int*)nonce)[i] = dist(gen);
    for (uint32_t i = 0; i < len % sizeof(int); ++i)
        nonce[i + len / sizeof(int)] = (uint8_t)dist(gen);
}

/**
 * @brief Add counter by len bytes
 *
 * @param ctr 16-byte couter array
 * @param len the number of bytes to be incremented
 */
void ctr128_inc_for_bytes(uint8_t* ctr, uint64_t len)
{
    if (ctr == NULL || len == 0)
        return;
    uint64_t inc_num = (len / 16) + (len % 16 ? 1 : 0);

    uint8_t rev_ctr[16];
    for (int i = 0; i < 16; ++i)
        rev_ctr[i] = ctr[15 - i];

    uint64_t* lo = (uint64_t*)(rev_ctr);
    uint64_t* hi = (uint64_t*)(rev_ctr + 8);

    uint64_t tmp_lo = *lo;
    *lo = *lo + inc_num;
    if (*lo < tmp_lo)
        *hi = *hi + 1;

    for (int i = 0; i < 16; ++i)
        ctr[i] = rev_ctr[15 - i];
}

/**
 * @brief The largest address aligned with 16 byte <= add for 128-bit AES CTR en/decryption
 */
inline uint32_t lowerBound(uint32_t keySize, uint32_t add)
{
    while (add * keySize & 15)
        --add;
    return add;
}

/**
 * @brief The smallest address aligned with 16 byte >=add,
 *  combined with lowerBound as the region to be en/decrypted.
 */
inline uint32_t upperBound(uint32_t keySize, uint32_t add)
{
    while (add * keySize & 15)
        ++add;
    return add;
}

static uint32_t OMP_BITONIC_CNT = 1 << 10;
/**
 * @brief Merge a bitonic sequence
 *
 * @tparam Key Key must support operator<()
 * @param a A decrypted array inside the enclave
 * @param low The starting index
 * @param cnt The number of elements
 * @param dir The direction, 1 is ascending, 0 is descending
 */
template <typename Key>
static void bitonicMerge(Key* arr, uint32_t low, uint32_t cnt, bool dir)
{
    if (cnt <= 1)
        return;
    if (arr == NULL)
        return;
    uint32_t k = 1;
    while (k < cnt)
        k <<= 1;
    k >>= 1;

    for (uint32_t i = low; i < low + cnt - k; ++i)
        if (dir == (arr[i] > arr[i + k]))
            std::swap(arr[i], arr[i + k]);

    bitonicMerge<Key>(arr, low, k, dir);
    bitonicMerge<Key>(arr, low + k, cnt - k, dir);
}

/**
 * @brief Merge an encrypted bitonic array with the help of secure buffer of the given size
 *
 * @tparam Key Key must support operator<()
 * @param secureBuf a buffer inside the enclave
 * @param secureBufSize secure buffer size
 * @param cArr The encrypted array to be merged
 * @param low The start index
 * @param cnt The number of elements
 * @param dir 1 ascending, 0 descending
 */
template <typename Key>
static void bitonicMergeEx(const sgx_aes_ctr_128bit_key_t* p_key,
    Key* secureBuf, size_t secureBufSize,
    Key* cArr, uint32_t low, uint32_t cnt, bool dir)
{
    if (cnt <= 1)
        return;
    if (secureBuf == NULL || cArr == NULL)
        return;
    uint32_t k = 1;
    while (k < cnt)
        k <<= 1;
    k >>= 1;

    uint32_t st = lowerBound(sizeof(Key), low);
    uint32_t ed = upperBound(sizeof(Key), low + cnt);

    // If the region fits into the secure buffer, transform into the previous one
    if ((ed - st) * sizeof(Key) <= secureBufSize)
    {
        uint8_t ctr_de[16] = { 0 };
        uint8_t ctr_en[16] = { 0 };
        ctr128_inc_for_bytes(ctr_de, sizeof(Key) * st);
        memcpy(ctr_en, ctr_de, 16);
        sgx_aes_ctr_decrypt(p_key, (uint8_t*)(cArr + st),
            (ed - st) * sizeof(Key), ctr_de, 128, (uint8_t*)secureBuf);

        low -= st;
#ifdef _OPENMP
#pragma omp parallel for if (cnt - k > OMP_BITONIC_CNT)
#endif
        for (uint32_t i = low; i < low + cnt - k; ++i)
            if (dir == (secureBuf[i] > secureBuf[i + k]))
                std::swap(secureBuf[i], secureBuf[i + k]);
#ifdef _OPENMP
#pragma omp taskgroup
        {
#pragma omp task  if (k > OMP_BITONIC_CNT)
#endif
            bitonicMerge<Key>(secureBuf, low, k, dir);
#ifdef _OPENMP
#pragma omp task  if (cnt - k > OMP_BITONIC_CNT)
#endif
            bitonicMerge<Key>(secureBuf, low + k, cnt - k, dir);
#ifdef _OPENMP
#pragma omp taskyield
        }
#endif
        sgx_aes_ctr_encrypt(p_key, (uint8_t*)secureBuf, (ed - st) * sizeof(Key),
            ctr_en, 128, (uint8_t*)(cArr + st));
        return;
    }

    // Enforce it to be multiple of 16
    const uint32_t batchCnt = lowerBound(sizeof(Key), secureBufSize / (2 * sizeof(Key)));
    ed = upperBound(sizeof(Key), low + cnt - k);
    uint32_t epoch = (ed - st) / batchCnt + ((ed - st) % batchCnt != 0);
    uint32_t rSt = lowerBound(sizeof(Key), low + k);
    uint32_t rEd = 0;
    Key* leftArr = cArr + st;
    Key* rightArr = cArr + rSt;
    uint32_t prevLeftLow = low;
    uint32_t prevRightLow = low + k;
    uint32_t _cnt = cnt;
    uint8_t ctr_lDe[16] = { 0 };
    uint8_t ctr_rDe[16] = { 0 };
    uint8_t ctr_lEn[16] = { 0 };
    uint8_t ctr_rEn[16] = { 0 };
    ctr128_inc_for_bytes(ctr_lDe, sizeof(Key) * st);
    ctr128_inc_for_bytes(ctr_rDe, sizeof(Key) * rSt);
    memcpy(ctr_lEn, ctr_lDe, 16);
    memcpy(ctr_rEn, ctr_rDe, 16);
    for (uint32_t e = 0; e < epoch; ++e)
    {
        uint32_t leftCnt = std::min(batchCnt, ed - st);
        leftCnt = lowerBound(sizeof(Key), leftCnt);
        uint32_t actualCnt = std::min(cnt - k, st + leftCnt - prevLeftLow);
        rEd = upperBound(sizeof(Key), prevRightLow + actualCnt);
        uint32_t rightCnt = rEd - rSt;

        sgx_aes_ctr_decrypt(p_key, (uint8_t*)leftArr,
            leftCnt * sizeof(Key), ctr_lDe, 128, (uint8_t*)secureBuf);
        sgx_aes_ctr_decrypt(p_key, (uint8_t*)rightArr,
            rightCnt * sizeof(Key), ctr_rDe, 128,
            (uint8_t*)(secureBuf + leftCnt));

        uint32_t lIdx = prevLeftLow - st;
        uint32_t rIdx = prevRightLow - rSt + leftCnt;
        // #ifdef _OPENMP
        // #pragma omp parallel for if (actualCnt > OMP_BITONIC_CNT)
        // #endif
        for (uint32_t i = 0; i < actualCnt; ++i)
            if (dir == (secureBuf[lIdx + i] > secureBuf[rIdx + i]))
                std::swap(secureBuf[lIdx + i], secureBuf[rIdx + i]);

        sgx_aes_ctr_encrypt(p_key, (uint8_t*)secureBuf,
            leftCnt * sizeof(Key), ctr_lEn, 128, (uint8_t*)leftArr);
        sgx_aes_ctr_encrypt(p_key, (uint8_t*)(secureBuf + leftCnt),
            rightCnt * sizeof(Key), ctr_rEn, 128, (uint8_t*)rightArr);

        leftArr += leftCnt;
        rightArr += rightCnt;
        st += leftCnt;
        rSt += rightCnt;
        prevLeftLow = st;
        prevRightLow = rSt;
        cnt -= actualCnt;
    }
    bitonicMergeEx<Key>(p_key, secureBuf, secureBufSize, cArr, low, k, dir);
    bitonicMergeEx<Key>(p_key, secureBuf, secureBufSize, cArr, low + k, _cnt - k, dir);
    return;
}

/**
 * @brief Sort a (small) decrypted array inside the enclave
 */
template <typename Key>
void bitonicSort(Key* arr, uint32_t low, uint32_t cnt, bool dir)
{
    if (cnt <= 1)
        return;
    if (arr == NULL)
        return;
    const uint32_t k = cnt >> 1;
#ifdef _OPENMP
#pragma omp taskgroup
    {
#pragma omp task  if (k > OMP_BITONIC_CNT)
#endif
        bitonicSort<Key>(arr, low, k, !dir);
#ifdef _OPENMP
#pragma omp task  if (cnt - k > OMP_BITONIC_CNT)
#endif
        bitonicSort<Key>(arr, low + k, cnt - k, dir);
#ifdef _OPENMP
#pragma omp taskyield
    }
#endif
    bitonicMerge<Key>(arr, low, cnt, dir);
}

/**
 * @brief Sort an encrypted array with the help of secure buffer of the given size
 *
 * @tparam Key Key must support operator<()
 * @param secureBuf a buffer inside the enclave
 * @param secureBufSize secure buffer size
 * @param cArr The encrypted array to be sorted
 * @param low The start index
 * @param cnt The number of elements
 * @param dir 1 ascending, 0 descending
 */
template <typename Key>
void bitonicSortEx(const sgx_aes_ctr_128bit_key_t* p_key,
    Key* secureBuf, size_t secureBufSize,
    Key* cArr, uint32_t low, uint32_t cnt, bool dir)
{
    if (cnt <= 1)
        return;
    if (secureBuf == NULL || cArr == NULL)
        return;
    const uint32_t k = cnt >> 1;
    const uint32_t st = lowerBound(sizeof(Key), low);
    const uint32_t ed = upperBound(sizeof(Key), low + cnt);
    if ((ed - st) * sizeof(Key) < secureBufSize)
    {
        uint8_t ctr[16] = { 0 };
        ctr128_inc_for_bytes(ctr, sizeof(Key) * st);
        sgx_aes_ctr_decrypt(
            p_key,
            (uint8_t*)(cArr + st), (ed - st) * sizeof(Key),
            ctr, 128,
            (uint8_t*)secureBuf
        );
        low -= st;
#ifdef _OPENMP
#pragma omp taskgroup
        {
#pragma omp task  if (k > OMP_BITONIC_CNT)
#endif
            bitonicSort<Key>(secureBuf, low, k, !dir);
#ifdef _OPENMP
#pragma omp task  if (cnt - k > OMP_BITONIC_CNT)
#endif
            bitonicSort<Key>(secureBuf, low + k, cnt - k, dir);
#ifdef _OPENMP
#pragma omp taskyield
        }
#endif
        bitonicMerge<Key>(secureBuf, low, cnt, dir);

        memset(ctr, 0, sizeof(ctr));
        ctr128_inc_for_bytes(ctr, sizeof(Key) * st);
        sgx_aes_ctr_encrypt(p_key, (uint8_t*)secureBuf, (ed - st) * sizeof(Key),
            ctr, 128, (uint8_t*)(cArr + st));
        return;
    }
    bitonicSortEx<Key>(p_key, secureBuf, secureBufSize, cArr, low, k, !dir);
    bitonicSortEx<Key>(p_key, secureBuf, secureBufSize, cArr, low + k, cnt - k, dir);
    bitonicMergeEx<Key>(p_key, secureBuf, secureBufSize, cArr, low, cnt, dir);
    return;
}

template <typename Key>
void OSort(const sgx_aes_ctr_128bit_key_t* p_key,
    Key* secureBuf, size_t secureBufSize,
    Key* cArr, uint32_t cnt, uint32_t threadNum)
{
    uint32_t bufferLoad = secureBufSize / sizeof(Key);
    OMP_BITONIC_CNT = std::min(cnt, bufferLoad) / threadNum;
#ifdef _OPENMP 
#pragma omp parallel num_threads(threadNum) if (cnt * sizeof(Key) * log2(cnt) > 6000000)
#pragma omp single
#endif
    bitonicSortEx<Key>(p_key, secureBuf, secureBufSize, cArr, 0, cnt, true);
}

// static std::atomic_int acnt; // for threading count

void bitonicMergeWithBuckets(
    uint8_t* secureBuf,
    uint32_t* shuffleKeys,
    uint32_t low, uint32_t cnt, bool dir,
    const sgx_aes_ctr_128bit_key_t* p_key,
    std::pair<uint32_t, uint32_t>* i_j,
    uint8_t* encryptedBuckets, uint32_t encBucketSize)
{
    if (cnt <= 1)
        return;
    uint32_t k = 1;
    while (k < cnt)
        k <<= 1;
    k >>= 1;
    static const uint32_t bucketSize = encBucketSize - SGX_AESGCM_IV_SIZE - SGX_AESGCM_MAC_SIZE;
    uint8_t* curLeftBucket = encryptedBuckets + (size_t)low * (size_t)encBucketSize;
    uint8_t* curRightBucket = encryptedBuckets + (size_t)(low + k) * (size_t)encBucketSize;
    // uint8_t* leftBucket = secureBuf + acnt * bucketSize * 2;
    // uint8_t* rightBucket = leftBucket + bucketSize;
    std::vector<uint8_t> leftBucketVec(bucketSize);
    std::vector<uint8_t> rightBucketVec(bucketSize);
    uint8_t* leftBucket = leftBucketVec.data();
    uint8_t* rightBucket = rightBucketVec.data();

    for (uint32_t i = low; i < low + cnt - k; ++i)
    {
        sgx_status_t ret = sgx_rijndael128GCM_decrypt(p_key,
            curLeftBucket + SGX_AESGCM_IV_SIZE, bucketSize,
            leftBucket,
            curLeftBucket, SGX_AESGCM_IV_SIZE,
            NULL, 0,
            (sgx_aes_gcm_128bit_tag_t*)(curLeftBucket + SGX_AESGCM_IV_SIZE + bucketSize)
        );
        if (ret != SGX_SUCCESS)
        {
            t_print("sgx_rijndael128GCM_decrypt left %d failed", i);
            return;
        }
        ret = sgx_rijndael128GCM_decrypt(p_key,
            curRightBucket + SGX_AESGCM_IV_SIZE, bucketSize,
            rightBucket,
            curRightBucket, SGX_AESGCM_IV_SIZE,
            NULL, 0,
            (sgx_aes_gcm_128bit_tag_t*)(curRightBucket + SGX_AESGCM_IV_SIZE + bucketSize)
        );
        if (ret != SGX_SUCCESS)
        {
            t_print("sgx_rijndael128GCM_decrypt right %d failed", i);
            return;
        }
        if (dir == (shuffleKeys[i] > shuffleKeys[i + k]))
        {
            std::swap(shuffleKeys[i], shuffleKeys[i + k]);
            std::swap(i_j[i], i_j[i + k]);
            std::swap(leftBucket, rightBucket);
        }

        // if (((Datum<int>*)leftBucket)[0].getDat() == 0 || ((Datum<int>*)rightBucket)[0].getDat() == 0)
        // {
        //     t_print("ERROR!!!! %d, %d\n", i, i + k);
        //     assert(false);
        // }

        genNonce(curLeftBucket, SGX_AESGCM_IV_SIZE);
        sgx_rijndael128GCM_encrypt(p_key,
            leftBucket, bucketSize,
            curLeftBucket + SGX_AESGCM_IV_SIZE,
            curLeftBucket, SGX_AESGCM_IV_SIZE,
            NULL, 0,
            (sgx_aes_gcm_128bit_tag_t*)(curLeftBucket + SGX_AESGCM_IV_SIZE + bucketSize)
        );
        if (ret != SGX_SUCCESS)
        {
            t_print("sgx_rijndael128GCM_encrypt left %d failed", i);
            return;
        }

        genNonce(curRightBucket, SGX_AESGCM_IV_SIZE);
        sgx_rijndael128GCM_encrypt(p_key,
            rightBucket, bucketSize,
            curRightBucket + SGX_AESGCM_IV_SIZE,
            curRightBucket, SGX_AESGCM_IV_SIZE,
            NULL, 0,
            (sgx_aes_gcm_128bit_tag_t*)(curRightBucket + SGX_AESGCM_IV_SIZE + bucketSize)
        );
        if (ret != SGX_SUCCESS)
        {
            t_print("sgx_rijndael128GCM_encrypt right %d failed", i);
            return;
        }
        curLeftBucket += encBucketSize;
        curRightBucket += encBucketSize;
    }
    uint32_t curCnt = 0;
#ifdef _OPENMP
#pragma omp taskgroup
    {
        // if (acnt < threadNum)
        // {
        //     ++curCnt;
        //     ++acnt;
#pragma omp task default(shared) if (cnt > OMP_BITONIC_CNT)
#endif
        bitonicMergeWithBuckets(secureBuf, shuffleKeys, low, k, dir,
            p_key, i_j, encryptedBuckets, encBucketSize);
#ifdef _OPENMP
        // }
        // else
        //     bitonicMergeWithBuckets(secureBuf, shuffleKeys, low, k, dir,
        //         p_key, i_j, encryptedBuckets, encBucketSize);

        // if (acnt < threadNum)
        // {
        //     ++curCnt;
        //     ++acnt;
#pragma omp task default(shared) if (cnt - k > OMP_BITONIC_CNT)
#endif
        bitonicMergeWithBuckets(secureBuf, shuffleKeys, low + k, cnt - k, dir,
            p_key, i_j, encryptedBuckets, encBucketSize);
#ifdef _OPENMP
        // }
        // else
        //     bitonicMergeWithBuckets(secureBuf, shuffleKeys, low + k, cnt - k, dir,
        //         p_key, i_j, encryptedBuckets, encBucketSize, threadNum);
#pragma omp taskyield
    }
#endif
    // acnt -= curCnt;
}

void bitonicSortWithBuckets(
    uint8_t* secureBuf,
    uint32_t* shuffleKeys,
    uint32_t low, uint32_t cnt, bool dir,
    const sgx_aes_ctr_128bit_key_t* p_key,
    std::pair<uint32_t, uint32_t>* i_j,
    uint8_t* encryptedBuckets, uint32_t encBucketSize)
{
    if (cnt <= 1)
        return;
    const uint32_t k = cnt >> 1;
    uint32_t curCnt = 0;

#ifdef _OPENMP
#pragma omp taskgroup
    {
        // if (acnt < threadNum)
        // {
        //     ++curCnt;
        //     ++acnt;
#pragma omp task default(shared) if (cnt > OMP_BITONIC_CNT)
#endif
        bitonicSortWithBuckets(secureBuf, shuffleKeys, low, k, !dir,
            p_key, i_j, encryptedBuckets, encBucketSize);

#ifdef _OPENMP
        // }
        // else
        //     bitonicSortWithBuckets(secureBuf, shuffleKeys, low, k, !dir,
        //         p_key, i_j, encryptedBuckets, encBucketSize);

        // if (acnt < threadNum)
        // {
        //     ++curCnt;
        //     ++acnt;
#pragma omp task default(shared) if (cnt - k > OMP_BITONIC_CNT)
#endif
        bitonicSortWithBuckets(secureBuf, shuffleKeys, low + k, cnt - k, dir,
            p_key, i_j, encryptedBuckets, encBucketSize);
#ifdef _OPENMP
        // }
        // else
        //     bitonicSortWithBuckets(secureBuf, shuffleKeys, low + k, cnt - k, dir,
        //         p_key, i_j, encryptedBuckets, encBucketSize, threadNum);
#pragma omp taskyield
    }
#endif

    // acnt -= curCnt;
    bitonicMergeWithBuckets(secureBuf, shuffleKeys, low, cnt, dir,
        p_key, i_j, encryptedBuckets, encBucketSize);
}

void bitonicSortWithBucketsParallel(
    uint8_t* secureBuf,
    uint32_t* shuffleKeys,
    uint32_t cnt,
    const sgx_aes_ctr_128bit_key_t* p_key,
    std::pair<uint32_t, uint32_t>* i_j,
    uint8_t* encryptedBuckets,
    uint32_t encBucketSize, uint32_t threadNum)
{
    OMP_BITONIC_CNT = cnt * 2 / threadNum + 1;
#ifdef _OPENMP
#pragma omp parallel num_threads(threadNum)
#pragma omp single
#endif
    bitonicSortWithBuckets(secureBuf, shuffleKeys, 0u, cnt, true,
        p_key, i_j, encryptedBuckets, encBucketSize);
}