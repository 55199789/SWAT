#ifndef DEFS_H
#define DEFS_H
#include <algorithm>
#include <array>
#include <cstddef>
#include <ctime>
#include <limits>
#include <string>

typedef unsigned char BYTE;  // 1byte
typedef unsigned short WORD; // 2bytes
typedef unsigned long DWORD; // 4bytes
#define MAX_PATH_LEN 4096
const size_t DEFAULT_SECURE_BUFFER_SIZE = 512 * 1024 * 1024;
#define RED "\033[0;31m"
#define GREEN "\033[0;32m"
#define YELLOW  "\033[0;33m"
#define BLUE  "\033[0;34m"
#define DEFAULT "\033[0m"
#define PAGE_SIZE sysconf(_SC_PAGESIZE)
#define DEFAULT_KEY_SIZE sizeof(uint32)
#define DEFAULT_PTR_SIZE sizeof(void *)
// #define OBLIVIOUSMERGE
// #define RECORD_LENGTH_IN_BYTE 8

#define DEFAULT_BUCKET_CAPACITY 256
#define DEFAULT_IN_PORT 12341
#define DEFAULT_THREAD_NUM 32
#define DEFAULT_EPSILON 1.0
#define DEFAULT_DELTA (1e-12)
#define DEFAULT_LAMBDA 512
#define DEFAULT_QUERY_NUM 100000
#define DEFAULT_K 8
#define DEFAULT_HOST "localhost"
// #define DEFAULT_OUT_PORT 3337
// #define DEFAULT_HOST_TYPE "rocksdb"
#define DEFAULT_OUT_PORT 6379
#define DEFAULT_HOST_TYPE "redis"
// #define DEFAULT_OUT_PORT 11211
// #define DEFAULT_HOST_TYPE "memcached"
#define DEFAULT_UPDATE_POLICY "increment"
#define DEFAULT_STORAGE_OVERHEAD 2.0
#define MIN_B 3
#define DEFAULT_STORAGE_CLIENT_NUM 32
#define DEFAULT_BATCH_SIZE 3
#define DEFAULT_NUM 0

#define CMD_LEN sizeof(OP)
#define REPLICA_DEL std::numeric_limits<uint32_t>::max()
enum class Types
{
    Int,
    LongLong,
    Float,
    Double,
    Record,
    TypeInvalid
};

Types resolveType(const std::string& type)
{
    if (type == "int")
        return Types::Int;
    if (type == "long long")
        return Types::LongLong;
    if (type == "float")
        return Types::Float;
    if (type == "double")
        return Types::Double;
    if (type == "record")
        return Types::Record;
    return Types::TypeInvalid;
}

enum class WeightUpdatePolicy
{
    Increment,
    Constant,
    Multiplication,
    UpdateTypeInvalid
};

WeightUpdatePolicy resolveUpdatePolicy(const std::string& type)
{
    if (type == "increment")
        return WeightUpdatePolicy::Increment;
    if (type == "constant")
        return WeightUpdatePolicy::Constant;
    if (type == "multiplication")
        return WeightUpdatePolicy::Multiplication;
    return WeightUpdatePolicy::UpdateTypeInvalid;
}

struct ServerParameter
{
    uint32_t inPort = DEFAULT_IN_PORT;
    uint32_t outPort = DEFAULT_OUT_PORT;
    uint32_t storageClientNum = DEFAULT_STORAGE_CLIENT_NUM;
    const char* storageHost = DEFAULT_HOST;
    const char* storageType = DEFAULT_HOST_TYPE;
};

struct ClientParameter
{
    uint32_t lambda = DEFAULT_LAMBDA;
    double epsilon = DEFAULT_EPSILON;
    double delta = DEFAULT_DELTA;
    uint32_t bucketCapacity = DEFAULT_BUCKET_CAPACITY;
    uint32_t k = DEFAULT_K;
    Types type = Types::Int;
    uint32_t seed = 11451;
    uint32_t threadNum = DEFAULT_THREAD_NUM;
    uint32_t batchSize = DEFAULT_BATCH_SIZE;
    uint32_t num = DEFAULT_NUM;
    uint32_t queryNum = DEFAULT_QUERY_NUM;
    double storageOverhead = DEFAULT_STORAGE_OVERHEAD;
    const char* serverPort = "12341";
    const char* serverName = DEFAULT_HOST;
    WeightUpdatePolicy updatePolicy = WeightUpdatePolicy::Constant;
    char* domainMin = "1";
    char* domainMax = "10000";
    double pointReads = 0.0;
    double smallReads = 0.0;
    double largeReads = 0.0;
    double inserts = 100.0;
    char* rangeLen = NULL;
};

template<typename Key>
struct BucketTag
{
    Key left, right;
};

std::string makeBucketID(uint8_t k, uint32_t i, uint32_t j = 0, uint64_t seed_ = 0)
{
    static uint64_t seed = seed_;
    static auto Cantor = [](uint32_t a, int32_t b) -> uint64_t
    {
        return uint64_t(a + b) * (a + b + 1ull) / 2 + b;
    };
    static auto mumurHash64 = [](uint64_t key, uint64_t seed) -> uint64_t
    {
        const uint64_t m = 0xc6a4a7935bd1e995;
        const int r = 47;

        uint64_t h = seed ^ m;

        uint64_t k = key;

        key *= m;
        key ^= key >> r;
        key *= m;

        h ^= key;
        h *= m;

        h ^= h >> r;
        h *= m;
        h ^= h >> r;

        return h;
    };

    static const char* digits = "0123456789abcdef";
    static const uint64_t lenDict[16] = {
        0x1ull,				  // 1
        0x10ull,			  // 2
        0x100ull,			  // 3
        0x1000ull,			  // 4
        0x10000ull,			  // 5
        0x100000ull,		  // 6
        0x1000000ull,		  // 7
        0x10000000ull,		  // 8
        0x100000000ull,		  // 9
        0x1000000000ull,	  // 10
        0x10000000000ull,	  // 11
        0x100000000000ull,	  // 12
        0x1000000000000ull,	  // 13
        0x10000000000000ull,  // 14
        0x100000000000000ull, // 15
        0x1000000000000000ull // 16
    };

    uint64_t hashVal = mumurHash64(Cantor(Cantor(k, i), j), seed);
    uint32_t len = std::upper_bound(lenDict, lenDict + 16, hashVal) - lenDict;

    std::string ret(len, '0');
    for (size_t i = 0, j = (len - 1) * 4; i < len; ++i, j -= 4)
        ret[i] = digits[(hashVal >> j) & 0x0f];

    return ret;
}

enum class OP
{
    SEARCH_OP
    , INSERT_OP
    , DELETE_OP
    // ,   PUSH_DOWN_OP
};

#pragma pack(16)
template <typename Key>
struct Datum
{
private:
    static int timestamp_()
    {
        static int timestamp_ = 0;
        return ++timestamp_;
    }

    Key dat;
    OP op;
    time_t timestamp;
    std::array<uint8_t, RECORD_LENGTH_IN_BYTE> val;

public:
    static Datum<Key> getDummy()
    {
        return Datum<Key>(std::numeric_limits<Key>::max(),
            OP::INSERT_OP, std::numeric_limits<time_t>::max());
    }

    Datum(const Key& dat_, const OP& op_,
        time_t ts = timestamp_(),
        std::array<uint8_t, RECORD_LENGTH_IN_BYTE> val_ = {})
        : dat(dat_),
        op(op_),
        timestamp(ts),
        val(val_)
    {
    }

    Datum() {}

    Datum(const Datum<Key>& another) : dat(another.dat),
        op(another.op),
        timestamp(another.timestamp),
        val(another.val)
    {
    }

    Datum(Datum<Key>&& another) : dat(std::exchange(another.dat, std::numeric_limits<Key>::max())),
        op(std::exchange(another.op, OP::INSERT_OP)),
        timestamp(std::exchange(another.timestamp, std::numeric_limits<time_t>::max())),
        val(std::move(another.val))
    {
    }

    void set(const Key& dat_, const OP& op_, time_t ts = timestamp_())
    {
        dat = dat_;
        op = op_;
        timestamp = ts;
    }

    const Key& getDat() const
    {
        return dat;
    }

    bool operator<(const Datum<Key>& another) const
    {
        if (dat < another.dat)
            return true;
        if (dat == another.dat)
            return timestamp > another.timestamp;
        return false;
    }

    bool operator>(const Datum<Key>& another) const
    {
        if (dat > another.dat)
            return true;
        if (dat == another.dat)
            return timestamp < another.timestamp;
        return false;
    }

    Datum<Key>& operator=(const Datum<Key>& another)
    {
        dat = another.dat;
        op = another.op;
        timestamp = another.timestamp;
        val = another.val;
        return *this;
    }

    inline size_t size() const
    {
        return sizeof(dat) + sizeof(val) + sizeof(op) + sizeof(timestamp);
    }
};
#pragma pack()

#endif // DEFS_H