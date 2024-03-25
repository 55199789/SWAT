#pragma once 

#include "libstorage/StorageInterface.h"
#include <assert.h>
#include <memory> 
#include <unordered_map>
#include <libmemcached/memcached.h>

class MemcachedStorage : public StorageInterface
{
public:
    MemcachedStorage() :
        client(memcached_create(NULL)) {}

    virtual ~MemcachedStorage() {}

    virtual void batchDel(const std::vector<std::string>& keys) override
    {
        if (keys.empty())
        {
            return;
        }

        for (const auto& key : keys)
        {
            memcached_return_t rc = memcached_delete(client.get(),
                key.c_str(), key.size(), 0);
            // delete keys from memcached 
            if (rc != MEMCACHED_SUCCESS)
            {
                throw std::runtime_error(std::string("Batch delete failed with ") +
                    key +
                    std::string(": ") +
                    memcached_strerror(client.get(), rc));
            }
        }

    }

    virtual void batchGet(const std::vector<std::string>& keys,
        uint8_t* dest) override
    {
        if (keys.empty())
        {
            return;
        }

        std::vector<const char*> cstrings;
        std::vector<size_t> keyLens;
        char return_key[MEMCACHED_MAX_KEY];
        size_t return_key_length;
        char* return_bucket;
        size_t return_bucket_length;
        uint32_t flags;

        // parse vector of strings to vector of cstrings and keyLens
        cstrings.reserve(keys.size());
        keyLens.reserve(keys.size());
        for (const auto& key : keys)
        {
            cstrings.push_back(key.c_str());
            keyLens.push_back(key.size());
        }

        // get keys from memcached 
        memcached_return_t rc = memcached_mget(client.get(),
            cstrings.data(), keyLens.data(), keys.size());
        if (rc != MEMCACHED_SUCCESS)
        {
            throw std::runtime_error(std::string("Batch read failed: ") +
                memcached_strerror(client.get(), rc));
        }

        while ((return_bucket = memcached_fetch(client.get(), return_key, &return_key_length,
            &return_bucket_length, &flags, &rc)))
        {
            if (rc != MEMCACHED_SUCCESS)
            {
                throw std::runtime_error(
                    std::string("Batch read failed inside while: ") +
                    memcached_strerror(client.get(), rc));
            }
            dest = std::copy(return_bucket, return_bucket + return_bucket_length, dest);
        }
    }

    virtual void batchSet(const std::vector<std::pair<std::string, std::string>>&
        keyValuePairs) override
    {
        if (keyValuePairs.empty())
        {
            return;
        }

        memcached_return rc;
        for (const auto& keyValue : keyValuePairs)
        {
            rc = memcached_set(client.get(),
                keyValue.first.c_str(), keyValue.first.size(),
                keyValue.second.c_str(), keyValue.second.size(),
                0, 0
            );
            if (rc != MEMCACHED_SUCCESS)
            {
                throw std::runtime_error(std::string("Batch set failed: ") + memcached_strerror(client.get(), rc));
            }
        }
    }

    virtual void setup(const char* ip, uint32_t port, uint32_t clientNum)
    {
        memcached_return rc;
        memcached_server_st* server = NULL;
        server = memcached_server_list_append(server, ip, port, &rc);
        rc = memcached_server_push(client.get(), server);
        if (rc != MEMCACHED_SUCCESS) {
            throw std::runtime_error(std::string("Batch read failed: ") +
                memcached_strerror(client.get(), rc));
        }
    }
private:
    std::unique_ptr<memcached_st> client;
};