#pragma once 

#include "libstorage/StorageInterface.h"
#include <omp.h>
#include <cpp_redis/cpp_redis>
#include <memory> 
#include <vector>

// class RedisStorage: public StorageInterface
// {
// public:
//     virtual ~RedisStorage() {}
//     virtual void batchDel(const std::vector<std::string>& keys) override
//     {
//         if (keys.size() < clients.size() * 2)
//         {
//             std::future<cpp_redis::reply> cmd = clients[0]->del(keys);
//             clients[0]->commit();
//             cpp_redis::reply reply = cmd.get();
//             if (reply.is_error())
//                 throw std::runtime_error("Batch del: " + reply.error());
//             return;
//         }

//         std::vector<std::future<cpp_redis::reply>> cmds;

//         uint32_t batchSize = keys.size() / clients.size()
//             + (0 != keys.size() % clients.size());

//         auto it = keys.begin();
//         for (uint32_t i = 0; i < clients.size(), it < keys.end(); ++i)
//         {
//             std::vector<std::string> curKey(it,
//                 std::min(keys.end(), it + batchSize));
//             cmds.push_back(
//                 clients[i]->del(curKey)
//             );
//             clients[i]->commit();

//             it += batchSize;
//         }

//         for (auto& cmd : cmds)
//         {
//             cpp_redis::reply reply = cmd.get();
//             if (reply.is_error())
//                 throw std::runtime_error("Batch del: " + reply.error());
//         }
//     }

//     virtual void batchGet(const std::vector<std::string>& keys,
//         uint8_t* dest) override
//     {
//         if (keys.size() < clients.size() * 2)
//         {
//             std::future<cpp_redis::reply> cmd = clients[0]->mget(keys);
//             clients[0]->commit();
//             cpp_redis::reply reply = cmd.get();
//             if (reply.is_error())
//                 throw std::runtime_error("Batch get: " + reply.error());

//             const std::vector<cpp_redis::reply>& replies = reply.as_array();
//             uint8_t* currentBucket = dest;
//             for (auto& reply : replies)
//             {
//                 const auto& bucketData = reply.as_string();
//                 currentBucket = std::copy(bucketData.begin(), bucketData.end(),
//                     currentBucket);
//             }
//             return;
//         }

//         std::vector<std::future<cpp_redis::reply>> cmds;

//         uint32_t batchSize = keys.size() / clients.size()
//             + (0 != keys.size() % clients.size());

//         auto it = keys.begin();
//         for (uint32_t i = 0; i < clients.size(), it < keys.end(); ++i)
//         {
//             std::vector<std::string> curKey(it,
//                 std::min(keys.end(), it + batchSize));
//             cmds.push_back(
//                 clients[i]->mget(curKey)
//             );
//             clients[i]->commit();

//             it += batchSize;
//         }

//         uint8_t* currentBucket = dest;
//         for (auto& cmd : cmds)
//         {
//             cpp_redis::reply reply = cmd.get();
//             if (reply.is_error())
//                 throw std::runtime_error("Batch get: " + reply.error());

//             const std::vector<cpp_redis::reply>& replies = reply.as_array();

//             for (auto& reply : replies)
//             {
//                 const auto& bucketData = reply.as_string();
//                 currentBucket = std::copy(bucketData.begin(), bucketData.end(),
//                     currentBucket);
//             }
//         }
//     }

//     virtual void batchSet(
//         const std::vector<std::pair<std::string, std::string>>& keyValuePairs) override
//     {
//         if (keyValuePairs.size() < clients.size() * 2)
//         {
//             std::future<cpp_redis::reply> cmd = clients[0]->mset(keyValuePairs);
//             clients[0]->commit();
//             cpp_redis::reply reply = cmd.get();
//             if (reply.is_error())
//                 throw std::runtime_error("Batch del: " + reply.error());
//             return;
//         }

//         std::vector<std::future<cpp_redis::reply>> cmds;

//         uint32_t batchSize = keyValuePairs.size() / clients.size()
//             + (0 != keyValuePairs.size() % clients.size());
//         auto it = keyValuePairs.begin();
//         for (uint32_t i = 0; i < clients.size(), it < keyValuePairs.end(); ++i)
//         {
//             std::vector<std::pair<std::string, std::string>> curKey(it,
//                 std::min(keyValuePairs.end(), it + batchSize));
//             cmds.push_back(
//                 clients[i]->mset(curKey)
//             );
//             clients[i]->commit();

//             it += batchSize;
//         }
//         for (auto& cmd : cmds)
//         {
//             cpp_redis::reply reply = cmd.get();
//             if (reply.is_error())
//                 throw std::runtime_error("Batch get: " + reply.error());
//         }
//     }

//     virtual void setup(const char* ip, uint32_t port, uint32_t clientNum) override
//     {
//         clients.clear();
//         clients.reserve(clientNum);
//         while (clientNum--)
//         {
//             clients.push_back(std::move(std::make_shared<cpp_redis::client>()));
//             clients.back()->connect(ip, port,
//                 [](const std::string& host,
//                     std::size_t port,
//                     cpp_redis::client::connect_state status) {
//                         if (status == cpp_redis::client::connect_state::dropped
//                             || status == cpp_redis::client::connect_state::failed
//                             || status == cpp_redis::client::connect_state::lookup_failed)
//                         {
//                             std::cerr
//                                 << "Redis client disconnected from "
//                                 << host << ":" << port
//                                 << std::endl;
//                             exit(-1);
//                         }
//                 });

//         }
//     }
// private:
//     std::vector<std::shared_ptr<cpp_redis::client>> clients;
// };

class RedisStorage : public StorageInterface
{
public:
    virtual ~RedisStorage() {}
    virtual void batchDel(const std::vector<std::string>& keys) override
    {
        if (keys.empty())
        {
            return;
        }

        std::future<cpp_redis::reply> cmd = client.del(keys);
        client.commit();
        cpp_redis::reply reply = cmd.get();
        if (reply.is_error())
            throw std::runtime_error("Batch del: " + reply.error());
        return;
    }

    virtual void batchGet(const std::vector<std::string>& keys,
        uint8_t* dest) override
    {
        if (keys.empty())
        {
            return;
        }

        std::future<cpp_redis::reply> cmd = client.mget(keys);
        client.commit();
        cpp_redis::reply reply = cmd.get();
        if (reply.is_error())
            throw std::runtime_error("Batch get: " + reply.error());

        const std::vector<cpp_redis::reply>& replies = reply.as_array();
        uint8_t* currentBucket = dest;
        for (auto& reply : replies)
        {
            const auto& bucketData = reply.as_string();
            currentBucket = std::copy(bucketData.begin(), bucketData.end(),
                currentBucket);
        }
        return;
    }

    virtual void batchSet(
        const std::vector<std::pair<std::string, std::string>>& keyValuePairs) override
    {
        if (keyValuePairs.empty())
        {
            return;
        }

        std::future<cpp_redis::reply> cmd = client.mset(keyValuePairs);
        client.commit();
        cpp_redis::reply reply = cmd.get();
        if (reply.is_error())
            throw std::runtime_error("Batch set: " + reply.error());
        return;
    }

    virtual void setup(const char* ip, uint32_t port, uint32_t clientNum) override
    {
        client.connect(ip, port,
            [](const std::string& host,
                std::size_t port,
                cpp_redis::client::connect_state status) {
                    if (status == cpp_redis::client::connect_state::dropped
                        || status == cpp_redis::client::connect_state::failed
                        || status == cpp_redis::client::connect_state::lookup_failed)
                    {
                        std::cerr
                            << "Redis client disconnected from "
                            << host << ":" << port
                            << std::endl;
                        exit(-1);
                    }
            });
    }
private:
    cpp_redis::client client;
};