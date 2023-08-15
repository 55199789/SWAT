#pragma once 

#include <string>
#include <vector>

enum class StorageType { REDIS, ROCKSDB, MEMCACHED, INVALID };

StorageType resolveStorageType(const std::string& type)
{
    if (type == "redis")
        return StorageType::REDIS;
    else if (type == "rocksdb")
        return StorageType::ROCKSDB;
    else if (type == "memcached")
        return StorageType::MEMCACHED;
    else
        return StorageType::INVALID;
}

class StorageInterface
{
public:
    virtual void batchDel(const std::vector<std::string>& keys) = 0;
    virtual void batchGet(const std::vector<std::string>& keys, uint8_t* dest) = 0;
    virtual void batchSet(const std::vector<std::pair<std::string, std::string>>&
        keyValuePairs) = 0;
    virtual void setup(const char* ip, uint32_t port, uint32_t clientNum) = 0;
    virtual ~StorageInterface() {}
};