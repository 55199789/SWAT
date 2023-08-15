#include "sgx_urts.h"
#include "sgx_uswitchless.h"
#include <assert.h>
#include <getopt.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <time.h>
#include "common/defs.h"
#include "libstorage/StorageInterface.h"
#include "libstorage/RedisStorage.hpp"
#include "libstorage/MemcachedStorage.hpp"
#include "libstorage/RocksStorage.hpp"
#include "server_u.h"

#include <iostream>

/* Global EID shared by multiple threads */
static sgx_enclave_id_t server_global_eid = 0;
static std::unique_ptr<StorageInterface> storageClient;

static std::vector<uint8_t> individualEncBucketsMemoryPool;
static std::vector<uint8_t> mergedBucketsMemoryPool;
static std::vector<uint8_t> encBucketsVec;

typedef struct _sgx_errlist_t
{
    sgx_status_t err;
    const char* msg;
    const char* sug; /* Suggestion */
} sgx_errlist_t;

/* Error code returned by sgx_create_enclave */
static sgx_errlist_t sgx_errlist[] = {
    {SGX_ERROR_UNEXPECTED,
     "Unexpected error occurred.",
     NULL},
    {SGX_ERROR_INVALID_PARAMETER,
     "Invalid parameter.",
     NULL},
    {SGX_ERROR_OUT_OF_MEMORY,
     "Out of memory.",
     NULL},
    {SGX_ERROR_ENCLAVE_LOST,
     "Power transition occurred.",
     "Please refer to the sample \"PowerTransition\" for details."},
    {SGX_ERROR_ENCLAVE_CRASHED,
    "SGX enclave crashed.",
     "Please refer to the sample \"EnclaveMemoryManagement\" for details."},
    {SGX_ERROR_INVALID_ENCLAVE,
     "Invalid enclave image.",
     NULL},
    {SGX_ERROR_INVALID_ENCLAVE_ID,
     "Invalid enclave identification.",
     NULL},
    {SGX_ERROR_INVALID_SIGNATURE,
     "Invalid enclave signature.",
     NULL},
    {SGX_ERROR_OUT_OF_EPC,
     "Out of EPC memory.",
     NULL},
    {SGX_ERROR_NO_DEVICE,
     "Invalid SGX device.",
     "Please make sure SGX module is enabled in the BIOS, and install SGX driver afterwards."},
    {SGX_ERROR_MEMORY_MAP_CONFLICT,
     "Memory map conflicted.",
     NULL},
    {SGX_ERROR_INVALID_METADATA,
     "Invalid enclave metadata.",
     NULL},
    {SGX_ERROR_DEVICE_BUSY,
     "SGX device was busy.",
     NULL},
    {SGX_ERROR_INVALID_VERSION,
     "Enclave version was invalid.",
     NULL},
    {SGX_ERROR_INVALID_ATTRIBUTE,
     "Enclave was not authorized.",
     NULL},
    {SGX_ERROR_ENCLAVE_FILE_ACCESS,
     "Can't open enclave file.",
     NULL},
};

/* Check error conditions for loading enclave */
void print_error_message(sgx_status_t ret)
{
    size_t idx = 0;
    size_t ttl = sizeof sgx_errlist / sizeof sgx_errlist[0];

    for (idx = 0; idx < ttl; idx++)
    {
        if (ret == sgx_errlist[idx].err)
        {
            if (NULL != sgx_errlist[idx].sug)
                printf("Info: %s\n", sgx_errlist[idx].sug);
            printf("Error: %s\n", sgx_errlist[idx].msg);
            break;
        }
    }

    if (idx == ttl)
        printf("Error code is 0x%X. Please refer to the \"Intel SGX SDK Developer Reference\" for more details.\n", ret);
}

sgx_status_t initialize_enclave()
{
    sgx_status_t resp = SGX_SUCCESS;
    sgx_uswitchless_config_t us_config = SGX_USWITCHLESS_CONFIG_INITIALIZER;
    us_config.num_uworkers = 9;
    us_config.num_tworkers = 0;
    const void* enclave_ex_p[32] = { 0 };

    enclave_ex_p[SGX_CREATE_ENCLAVE_EX_SWITCHLESS_BIT_IDX] = (const void*)&us_config;
    resp = sgx_create_enclave_ex(ENCLAVE_PATH, SGX_DEBUG_FLAG, NULL, NULL, &server_global_eid,
        NULL, SGX_CREATE_ENCLAVE_EX_SWITCHLESS, enclave_ex_p);

    if (resp != SGX_SUCCESS)
    {
        print_error_message(resp);
        return resp;
    }

    return resp;
}

void ocall_setSeed(uint32_t seed)
{
    makeBucketID(0, 0, 0, seed);
}

void ocall_timing(uint8_t start, const char* msg)
{
    static auto st = std::chrono::high_resolution_clock::now();
    static auto ed = std::chrono::high_resolution_clock::now();
    if (start == true)
    {
        st = std::chrono::high_resolution_clock::now();
    }
    else
    {
        ed = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(ed - st);
        printf("%s time: %ld microseconds\n", msg, duration.count());
    }
}

void terminate_enclave()
{
    sgx_status_t result = sgx_destroy_enclave(server_global_eid);
    if (result != SGX_SUCCESS)
    {
        printf("Failed to terminate enclave.\n");
    }
    else
    {
        printf("Enclave successfully terminated.\n");
    }
}


void ocall_resizeMemoryPool(uint8_t** ptr, uint32_t poolIndex, size_t size)
{
    switch (poolIndex)
    {
    case 0:
        if (size > individualEncBucketsMemoryPool.size())
            individualEncBucketsMemoryPool.resize(size);
        *ptr = individualEncBucketsMemoryPool.data();
        break;
    case 1:
        if (size > mergedBucketsMemoryPool.size())
            mergedBucketsMemoryPool.resize(size);
        *ptr = mergedBucketsMemoryPool.data();
        break;
    default:
        *ptr = NULL;
        printf("Invalid pool index %d\n \
                \t0: individualEncBucketsMemoryPool\n \
                \t1: mergedBucketsMemoryPool\n", poolIndex);
        break;
    }
}

void ocall_getBucketData(const uint32_t* replicaCnts, size_t rCnt,
    const uint32_t* bucketCnts, size_t bCnt,
    uint8_t** encBucket,
    uint32_t encBucketSize)
{
    uint32_t totalBuckets = 0;
    std::vector<std::string> bucketGetIDStr;
    std::vector<std::string> bucketDelIDStr;
    for (uint32_t lvl = 1; lvl < bCnt; ++lvl)
    {
        for (uint32_t i = 0; i < bucketCnts[lvl]; ++i)
        {
            bucketGetIDStr.push_back(makeBucketID(lvl - 1, i));
            for (uint32_t j = 0; j < replicaCnts[i + totalBuckets]; ++j)
            {
                bucketDelIDStr.push_back(makeBucketID(lvl - 1, i, j));
            }
        }
        for (uint32_t j = 0; j < replicaCnts[bucketCnts[lvl] + totalBuckets]; ++j)
        {
            bucketDelIDStr.push_back(makeBucketID(lvl - 1, bucketCnts[lvl], j));
        }
        totalBuckets += bucketCnts[lvl] + 1;
    }

    totalBuckets += bucketCnts[0];
    mergedBucketsMemoryPool.resize(totalBuckets * encBucketSize);
    storageClient->batchGet(bucketGetIDStr,
        mergedBucketsMemoryPool.data() + (encBucketSize - 28) * bucketCnts[0]);
    *encBucket = mergedBucketsMemoryPool.data();
    storageClient->batchDel(bucketDelIDStr);
}

void ocall_insertBuckets(uint8_t lvl, uint32_t encBucketSize,
    const uint8_t* individualEncBuckets,
    const uint8_t* i_j_, size_t size_)
{
    uint32_t cnt = size_ / sizeof(std::pair<uint32_t, uint32_t>);
    const uint32_t batchSize = std::min(cnt, (500 << 20) / encBucketSize);
    std::vector<std::pair<std::string, std::string>> bucketIDDataPairs(batchSize);
    std::pair<uint32_t, uint32_t>* i_j = (std::pair<uint32_t, uint32_t>*)i_j_;
    const uint8_t* curBucket = individualEncBuckets;
    for (uint32_t i = 0; i < cnt; ++i)
    {
        std::string encBucket(curBucket, curBucket + encBucketSize);
        bucketIDDataPairs[i % batchSize] = std::make_pair(
            makeBucketID(lvl, i_j[i].first, i_j[i].second),
            encBucket
        );
        if ((i + 1) % batchSize == 0)
        {
            storageClient->batchSet(bucketIDDataPairs);
        }
        assert(curBucket < individualEncBucketsMemoryPool.data() + individualEncBucketsMemoryPool.size());
        curBucket += encBucketSize;
    }
    if (cnt % batchSize != 0)
    {
        bucketIDDataPairs.resize(cnt % batchSize);
        storageClient->batchSet(bucketIDDataPairs);
    }
}

void ocall_getBuckets(uint8_t* bucketIDs_, size_t size_, uint8_t** encBuckets, uint32_t encBucketSize, uint32_t* bucketCnt)
{
    if (bucketIDs_ == NULL)
    {
        printf("ocall_getBuckets bucketIDs is NULL\n");
        return;
    }
    std::vector<std::string> bucketIDs;
    char* curBucketID = strtok((char*)(bucketIDs_), " ");
    while (curBucketID != NULL)
    {
        ++(*bucketCnt);
        bucketIDs.emplace_back(curBucketID);
        curBucketID = strtok(NULL, " ");
    }
    // get bucket data from redis
    if (size_t(*bucketCnt) * encBucketSize > encBucketsVec.size())
    {
        encBucketsVec.resize((size_t)(*bucketCnt) * (size_t)encBucketSize);
    }
    *encBuckets = encBucketsVec.data();
    // *encBuckets = new uint8_t[(*bucketCnt) * encBucketSize];
    storageClient->batchGet(bucketIDs, *encBuckets);
}

void ocall_delete(void* addr)
{
    // individualEncBucketsMemoryPool.clear();
    // mergedBucketsMemoryPool.clear();
    // encBucketsVec.clear();
    // delete[](uint8_t*)addr;
}

void usage()
{
    std::cout << "Server proxy: \n";
    // Network Parameters
    std::cout << "\t --storageHost: Storage server host name\n";
    std::cout << "\t --storageType: Storage server type (redis, memcached, rocksdb)\n";
    std::cout << "\t --outPort: Storage server port\n";
    std::cout << "\t --inPort: Input port\n";
    std::cout << "\t --storageClientNum: Number of clients for the KV store\n";
};

int main(int argc, char* argv[])
{
    sgx_status_t result = SGX_SUCCESS;
    int ret = 1;
    ServerParameter par;

    int o;
    int option_index = 0;
    struct option longopts[] = {
        {"storageHost", optional_argument, NULL, 'h'},
        {"storageType", optional_argument, NULL, 't'},
        {"inPort", optional_argument, NULL, 'i'},
        {"outPort", optional_argument, NULL, 'o'},
        {"storageClientNum", optional_argument, NULL, 's'},
        {NULL, 0, NULL, 0} };

    while ((o = getopt_long(argc, argv, "hio",
        longopts, &option_index)) != -1)
    {
        switch (o)
        {
        case 'h':
            par.storageHost = optarg;
            break;
        case 't':
            par.storageType = optarg;
            break;
        case 'i':
            par.inPort = std::atoi(optarg);
            break;
        case 'o':
            par.outPort = std::atoi(optarg);
            break;
        case 'c':
            par.storageClientNum = std::atoi(optarg);
        default:
            usage();
            exit(-1);
        }
    }

    printf("Host: Creating an tls server enclave\n");
    result = initialize_enclave();
    if (result != SGX_SUCCESS)
    {
        goto exit;
    }
    printf("Host: The tls server enclave has been created\n");

    switch (resolveStorageType(par.storageType))
    {
    case StorageType::REDIS:
        storageClient = std::unique_ptr<RedisStorage>(new RedisStorage());
        break;
    case StorageType::ROCKSDB:
        storageClient = std::unique_ptr<RocksStorage>(new RocksStorage());
        break;
    case StorageType::MEMCACHED:
        storageClient = std::unique_ptr<MemcachedStorage>(new MemcachedStorage());
        break;
    default:
        throw std::runtime_error("unknown storage type");
        break;
    }

    storageClient->setup(par.storageHost, par.outPort, par.storageClientNum);

    printf("calling setup_server\n");
    result = set_up_server(server_global_eid, &ret, &par);
    if (result != SGX_SUCCESS || ret != 0)
    {
        printf("Host: setup_server failed\n");
        print_error_message(result);
        goto exit;
    }

exit:

    terminate_enclave();
    return ret;
}
