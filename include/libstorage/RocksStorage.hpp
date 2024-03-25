#pragma once
#include "libstorage/StorageInterface.h"
#include <memory>
#include <rocksdb/db.h>

class RocksStorage : public StorageInterface
{
public:
    RocksStorage()
    {
        rocksdb::Options options;
        // Optimize RocksDB. This is the easiest way to get RocksDB to perform well
        options.IncreaseParallelism();
        options.OptimizeLevelStyleCompaction();
        // create the DB if it's not already present
        options.create_if_missing = true;
        rocksdb::Status status =
            rocksdb::DB::Open(options, "/tmp/swatdb", &db);
        if (!status.ok())
        {
            throw std::runtime_error("RocksDB open failed");
        }
    }

    virtual ~RocksStorage()
    {
        delete db;
    }

    virtual void batchDel(const std::vector<std::string> &keys) override
    {
        if (keys.empty())
        {
            return;
        }

        rocksdb::WriteBatch batch;
        for (const auto &key : keys)
        {
            batch.Delete(key);
        }
        rocksdb::Status status = db->Write(wOpt, &batch);
        if (!status.ok())
        {
            throw std::runtime_error("RocksDB batch delete failed");
        }
    }

    virtual void batchGet(const std::vector<std::string> &keys, uint8_t *dest) override
    {
        if (keys.empty())
        {
            return;
        }

        for (const std::string &key : keys)
        {
            std::string value;
            rocksdb::Status status = db->Get(rOpt, key, &value);
            if (!status.ok())
            {
                throw std::runtime_error("RocksDB get failed");
            }
            dest = std::copy(value.begin(), value.end(), dest);
        }
    }
    virtual void batchSet(const std::vector<std::pair<std::string, std::string>> &
                              keyValuePairs) override
    {
        if (keyValuePairs.empty())
        {
            return;
        }

        rocksdb::WriteBatch batch;
        for (const auto &kv : keyValuePairs)
        {
            batch.Put(kv.first, kv.second);
        }
        rocksdb::Status status = db->Write(wOpt, &batch);
        if (!status.ok())
        {
            throw std::runtime_error("RocksDB write failed");
        }
    }
    virtual void setup(const char *ip, uint32_t port, uint32_t clientNum) override
    {
    }

private:
    rocksdb::DB *db;
    rocksdb::WriteOptions wOpt;
    rocksdb::ReadOptions rOpt;
};