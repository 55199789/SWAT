/**
 *
 * MIT License
 *
 * Copyright (c) Open Enclave SDK contributors.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE
 *
 */

#include <arpa/inet.h>
#include <errno.h>
#include <getopt.h>
#include <netdb.h>
#include <netinet/in.h>
#include <resolv.h>
#include <sys/socket.h>
#include <string.h>
#include <unistd.h>

#include <atomic>
#include <chrono>
#include <iostream>
#include <memory>
#include <mutex>
#include <sstream>
#include <string>
#include <thread>
#include <unordered_map>
#include <type_traits>
#include <vector>

#include "sgx_utls.h"

#include <openssl/bio.h>
#include <openssl/err.h>
#include <openssl/pem.h>
#include <openssl/ssl.h>
#include <openssl/x509.h>
#include <openssl/x509_vfy.h>

#include "common/common.h"
#include "common/defs.h"
#include "Query.hpp"
#include "UniformQueryFactory.hpp"
#include "SamplePool.hpp"

int verify_callback(int preverify_ok, X509_STORE_CTX *ctx);

void usage()
{
    std::cout << "Client proxy: \n"
              << "\t --serverName: Server host name\n"
              << "\t --serverPort: Server port\n"
              << "\t --lambda: Security parameter\n"
              << "\t --epsilon: Differential obliviousness budget\n"
              << "\t --delta: Failure probability\n"
              << "\t --seed: Random seed for reproducing\n"
              << "\t --k: Number of layers for dynamization\n"
              << "\t --type: Type of keys\n"
              << "\t --storageOverhead: Storage overhead\n"
              << "\t --batch: Batch size for frequency smoothign\n"
              << "\t --updatePolicy: Update policy for frequency smoothing\n"
              << "\t --domainMin: Minimum value of the domain\n"
              << "\t --domainMax: Maximum value of the domain\n"
              << "\t --rangeLen: Maximum length of the range\n"
              << "\t --pointReads: Number of point reads\n"
              << "\t --smallReads: Number of small reads\n"
              << "\t --largeReads: Number of large reads\n"
              << "\t --inserts: Number of inserts\n"
              << "\t --num: Number of initial data\n"
              << "\t --threadNum: Number of threads\n"
              << std::endl;
};

template <typename T>
std::ostream &operator<<(typename std::enable_if<std::is_enum<T>::value, std::ostream>::type &stream, const T &e)
{
    return stream << static_cast<typename std::underlying_type<T>::type>(e);
}

int parse_arguments(
    int argc,
    char **argv,
    ClientParameter &par)
{
    int o;
    int option_index = 0;
    struct option longopts[] = {
        {"serverName", optional_argument, NULL, 'h'},
        {"serverPort", optional_argument, NULL, 'p'},
        {"lambda", optional_argument, NULL, 'l'},
        {"epsilon", optional_argument, NULL, 'e'},
        {"bucketCapacity", optional_argument, NULL, 'B'},
        {"seed", optional_argument, NULL, 's'},
        {"k", optional_argument, NULL, 'k'},
        {"type", optional_argument, NULL, 't'},
        {"storageOverhead", optional_argument, NULL, 'o'},
        {"batch", optional_argument, NULL, 'b'},
        {"updatePolicy", optional_argument, NULL, 'u'},
        {"domainMin", optional_argument, NULL, 'm'},
        {"domainMax", optional_argument, NULL, 'M'},
        {"pointReads", optional_argument, NULL, 'P'},
        {"smallReads", optional_argument, NULL, 'r'},
        {"largeReads", optional_argument, NULL, 'R'},
        {"rangeLen", optional_argument, NULL, 'L'},
        {"inserts", optional_argument, NULL, 'i'},
        {"num", optional_argument, NULL, 'n'},
        {"queryNum", optional_argument, NULL, 'N'},
        {"threadNum", optional_argument, NULL, 'f'},
        {NULL, 0, NULL, 0}};

    while ((o = getopt_long(argc, argv, "hpledsktobumML",
                            longopts, &option_index)) != -1)
    {
        switch (o)
        {
        case 'h':
            par.serverName = optarg;
            break;
        case 'p':
            par.serverPort = optarg;
            break;
        case 'l':
            par.lambda = std::atoi(optarg);
            break;
        case 'e':
            par.epsilon = std::atof(optarg);
            break;
        case 'B':
            par.bucketCapacity = std::atoi(optarg);
            break;
        case 'd':
            par.delta = std::atoi(optarg);
            break;
        case 's':
            par.seed = std::atoi(optarg);
            break;
        case 'k':
            par.k = std::atoi(optarg);
            break;
        case 't':
            par.type = resolveType(std::string(optarg));
            break;
        case 'o':
            par.storageOverhead = std::atof(optarg);
            // assert(par.storageOverhead >= 2.0);
            break;
        case 'b':
            par.batchSize = std::atoi(optarg);
            break;
        case 'u':
            par.updatePolicy = resolveUpdatePolicy(std::string(optarg));
            break;
        case 'm':
            par.domainMin = optarg;
            break;
        case 'M':
            par.domainMax = optarg;
            break;
        case 'P':
            par.pointReads = std::atof(optarg);
            break;
        case 'r':
            par.smallReads = std::atof(optarg);
            break;
        case 'R':
            par.largeReads = std::atof(optarg);
            break;
        case 'i':
            par.inserts = std::atof(optarg);
            break;
        case 'L':
            par.rangeLen = optarg;
            break;
        case 'n':
            par.num = std::atoi(optarg);
            break;
        case 'N':
            par.queryNum = std::atoi(optarg);
            break;
        case 'f':
            par.threadNum = std::atoi(optarg);
            break;
        default:
            usage();
            return -1;
        }
    }

    par.num = par.num / par.bucketCapacity * par.bucketCapacity;
    par.delta = exp(-log(par.lambda) * log(par.lambda));
    std::cout
        << "lambda: " << par.lambda
        << ", epsilon: " << par.epsilon
        << ", bucketCapacity: " << par.bucketCapacity
        << ", seed: " << par.seed
        << ", k: " << par.k
        << ", type: " << par.type
        << ", storageOverhead: " << par.storageOverhead
        << ", batchSize: " << par.batchSize
        << ", updatePolicy: " << par.updatePolicy
        << ", domainMin: " << par.domainMin
        << ", domainMax: " << par.domainMax
        << ", pointReads: " << par.pointReads
        << ", smallReads: " << par.smallReads
        << ", largeReads: " << par.largeReads
        << ", rangeLen: " << (par.rangeLen == NULL ? "" : par.rangeLen)
        << ", inserts: " << par.inserts
        << ", num: " << par.num
        << ", queryNum: " << par.queryNum
        << ", threadNum: " << par.threadNum
        << ", recordLength: " << RECORD_LENGTH_IN_BYTE
        << std::endl;
    return 0;
}

// create a socket and connect to the server_name:server_port
int create_socket(const char *server_name, const char *server_port)
{
    int sockfd = -1;
    struct addrinfo hints, *dest_info, *curr_di;
    int res;

    hints = {0};
    hints.ai_family = AF_INET;
    hints.ai_socktype = SOCK_STREAM;

    if ((res = getaddrinfo(server_name, server_port, &hints, &dest_info)) != 0)
    {
        printf(
            CLIENT_PROXY "Error: Cannot resolve hostname %s. %s\n",
            server_name,
            gai_strerror(res));
        goto done;
    }

    curr_di = dest_info;
    while (curr_di)
    {
        if (curr_di->ai_family == AF_INET)
        {
            break;
        }

        curr_di = curr_di->ai_next;
    }

    if (!curr_di)
    {
        printf(
            CLIENT_PROXY "Error: Cannot get address for hostname %s.\n",
            server_name);
        goto done;
    }

    sockfd = socket(AF_INET, SOCK_STREAM, 0);
    if (sockfd == -1)
    {
        printf(CLIENT_PROXY "Error: Cannot create socket %d.\n", errno);
        goto done;
    }

    if (connect(
            sockfd,
            (struct sockaddr *)curr_di->ai_addr,
            sizeof(struct sockaddr)) == -1)
    {
        printf(
            CLIENT_PROXY "failed to connect to %s:%s (errno=%d)\n",
            server_name,
            server_port,
            errno);
        close(sockfd);
        sockfd = -1;
        goto done;
    }
    printf(CLIENT_PROXY "connected to %s:%s\n", server_name, server_port);

done:
    if (dest_info)
        freeaddrinfo(dest_info);

    return sockfd;
}

size_t readFromServer(SSL *ssl, uint8_t *buffer, size_t bufferSize)
{
    size_t accBytesRead = 0;
    do
    {
        uint32_t bytesRead = SSL_read(ssl, buffer, bufferSize);

        if (bytesRead <= 0)
        {
            int error = SSL_get_error(ssl, bytesRead);
            if (error == SSL_ERROR_WANT_READ)
                continue;

            return accBytesRead;
        }

        accBytesRead += bytesRead;
        buffer += bytesRead;
    } while (accBytesRead < bufferSize);
    return accBytesRead;
}

void writeToServer(SSL *ssl, uint8_t *data, size_t size)
{
    int error = 0;
    int bytes_written = 0;

    while ((bytes_written = SSL_write(ssl, data, size)) <= 0)
    {
        error = SSL_get_error(ssl, bytes_written);
        if (error == SSL_ERROR_WANT_WRITE)
            continue;
        printf(CLIENT_PROXY "Failed to write buffer! SSL_write returned %d\n", error);
        printf(CLIENT_PROXY "error: %s\n", strerror(errno));
        exit(-1);
    }
}

template <typename T>
void uploadToServer(SSL *ssl, const std::vector<T> &buffer)
{
    const size_t bufferSize = buffer.size() * sizeof(T);
    writeToServer(ssl, (uint8_t *)&bufferSize, sizeof(bufferSize));
    writeToServer(ssl, (uint8_t *)buffer.data(), bufferSize);
}

// This routine conducts a simple HTTP request/response communication with
// server
int sendParameters(SSL *ssl, ClientParameter &par)
{
    unsigned char buf[50];
    int ret = 1;
    int error = 0;
    unsigned int len = 0;
    int bytes_written = 0;
    int bytes_read = 0;
    uint32_t binCapacity = 0;

    printf(CLIENT_PROXY "write parameters to the server:\n");
    while ((bytes_written = SSL_write(ssl, &par, sizeof(par))) <= 0)
    {
        error = SSL_get_error(ssl, bytes_written);
        if (error == SSL_ERROR_WANT_WRITE)
            continue;
        printf(CLIENT_PROXY "Failed! SSL_write returned %d\n", error);
        ret = bytes_written;
        goto done;
    }
    printf(CLIENT_PROXY "%d bytes written\n", bytes_written);

    readFromServer(ssl, (uint8_t *)&binCapacity, sizeof(binCapacity));
    printf(CLIENT_PROXY "binCapacity: %d\n", binCapacity);

    ret = 0;
done:
    return ret;
}

template <typename Key>
void updateSamplePool(SSL *ssl)
{
    uint64_t sizeLevel = 0;
    uint32_t size;
    uint32_t level;
    int error = 0;
    int bytes_written = 0;

    readFromServer(ssl, (uint8_t *)&sizeLevel, sizeof(sizeLevel));
    size = sizeLevel >> 32;
    level = sizeLevel & 0xFFFFFFFF;
    std::vector<BucketTag<Key>> bucketTagList(size / sizeof(BucketTag<Key>));

    readFromServer(ssl, (uint8_t *)bucketTagList.data(), size);
    std::vector<uint32_t> replicaCnt = SamplePool<Key>::getInstance().update(level, bucketTagList);
    uploadToServer(ssl, replicaCnt);
}

template <typename Key>
int generateQueries(SSL *ssl, ClientParameter &par)
{
    // we use int here
    const size_t IST_BUFFER_LEN = par.bucketCapacity * sizeof(Datum<Key>);
    // const size_t SETUP_BUFFER_LEN = CMD_LEN + par.bucketCapacity * sizeof(Datum<Key>);
    std::mutex opMutex;
    std::vector<uint8_t> insertBufferCMD(IST_BUFFER_LEN);
    Datum<Key> *insertBuffer = (Datum<Key> *)(insertBufferCMD.data());
    // Datum<Key>* setupBuffer = (Datum<Key>*) (setupBufferCMD.data() + CMD_LEN );
    Key min;
    Key max;
    bool pureRangeQuery = false;
    Key rangeLen;
    std::istringstream(par.domainMin) >> min;
    std::istringstream(par.domainMax) >> max;
    printf("rangeLen: %s\n", par.rangeLen);
    if (par.rangeLen != NULL)
    {
        std::istringstream(par.rangeLen) >> rangeLen;
        pureRangeQuery = true;
    }
    else
        rangeLen = max - min;
    std::cout << "Debug\tmin: " << min << ", max: " << max << ", rangeLen: " << rangeLen << std::endl;
    SamplePool<Key>::getInstance().setup(
        par.k + 1,
        par.batchSize,
        par.seed,
        par.storageOverhead,
        min, max, rangeLen,
        par.updatePolicy);
    uint32_t bufferCnt = 0;

    int ret = 1;
    time_t cnt = 0;
    OP insertOP = OP::INSERT_OP;
    // OP pushDownOP = OP::PUSH_DOWN_OP;
    auto st = std::chrono::high_resolution_clock::now();
    auto ed = std::chrono::high_resolution_clock::now();
    Key smallRangeLen = 100.0 / (par.num + 1) * (max - min);
    Key largeRangLen = 0.05 * (max - min);
    /***
     *  Setup in one phase
     */
    if (par.num > 0)
    {
        st = std::chrono::high_resolution_clock::now();
        size_t initDBSize = par.num * sizeof(Datum<Key>);
        writeToServer(ssl, (uint8_t *)&insertOP, sizeof(insertOP));
        writeToServer(ssl, (uint8_t *)&initDBSize, sizeof(initDBSize));
        // memcpy(setupBufferCMD.data() + sizeof(OP), &par.num, sizeof(uint32_t));
        UniformQueryFactory<Key>::getInstance().setup(par.seed, min, max,
                                                      smallRangeLen, largeRangLen, {0.0, 0.0, 0.0, 1.0, 0.0});
        for (uint32_t i = 0; i < par.num; ++i)
        {
            std::unique_ptr<Query<Key>> curQuery = UniformQueryFactory<Key>::getInstance().createQuery();
            switch (curQuery->getOp())
            {
            case OP::INSERT_OP:
            case OP::DELETE_OP:
                insertBuffer[bufferCnt++].set(curQuery->getQueryData().first, curQuery->getOp(), ++cnt);
                break;
            case OP::SEARCH_OP:
            default:
                throw std::runtime_error("unknown operation");
                break;
            }
            if (bufferCnt == par.bucketCapacity)
            {
                writeToServer(ssl, (uint8_t *)insertBuffer, IST_BUFFER_LEN);
                bufferCnt = 0;
            }
        }
        if (bufferCnt > 0)
        {
            writeToServer(ssl, (uint8_t *)insertBuffer, bufferCnt * sizeof(Datum<Key>));
            bufferCnt = 0;
        }
        updateSamplePool<Key>(ssl);
        readFromServer(ssl, (uint8_t *)&ret, sizeof(ret));
        if (ret != 0)
        {
            printf(CLIENT_PROXY "Error in setup phase!\n");
            exit(-1);
        }
        ed = std::chrono::high_resolution_clock::now();
        printf("Setting up %d keys took %" PRId64 " milliseconds\n",
               par.num,
               std::chrono::duration_cast<std::chrono::milliseconds>(ed - st).count());
    }
    // setupBufferCMD.clear();
    // printf("Input Enter to Continue...\n");
    // getchar();

    /***
     *  Inserting one by one
     */
    // memcpy(insertBufferCMD.data(), &insertOP, sizeof(OP));
    // bufferCnt = 0;
    // UniformQueryFactory<Key>::getInstance().setup(par.seed, min, max,
    // smallRangeLen, largeRangLen, { 0.0, 0.0, 0.0, 1.0, 0.0 });
    // st = std::chrono::high_resolution_clock::now();
    // for (uint32_t i = 0; i < par.num; ++i)
    // {
    //     std::unique_ptr<Query<Key>> curQuery = UniformQueryFactory<Key>::getInstance().createQuery();
    //     switch (curQuery->getOp())
    //     {
    //     case OP::INSERT_OP:
    //     case OP::DELETE_OP:
    //         insertBuffer[bufferCnt++].set(curQuery->getQueryData().first, curQuery->getOp(), ++cnt);
    //         break;
    //     case OP::SEARCH_OP:
    //     default:
    //         throw std::runtime_error("unknown operation");
    //         break;
    //     }

    //     if (bufferCnt == par.bucketCapacity)
    //     {
    //         bufferCnt = 0;
    //         uploadToServer(ssl, insertBufferCMD);
    //         updateSamplePool<Key>(ssl);
    //     }
    // }
    // ed = std::chrono::high_resolution_clock::now();
    // printf("Inserting %d keys took %" PRId64 " milliseconds\n",
    //     par.num, std::chrono::duration_cast<std::chrono::milliseconds>(ed - st).count());
    // printf("Input Enter to Continue...\n");
    // getchar();

    /***
     * PUSH DOWN
     */
    // st = std::chrono::high_resolution_clock::now();
    // memcpy(insertBufferCMD.data(), &pushDownOP, sizeof(OP));
    // uploadToServer(ssl, insertBufferCMD);
    // updateSamplePool<Key>(ssl);
    // readFromServer(ssl, &ret, sizeof(ret));
    // if (ret != 0)
    // {
    //     printf(CLIENT_PROXY "Error in setup phase!\n");
    //     exit(-1);
    // }
    // ed = std::chrono::high_resolution_clock::now();
    // printf("Pushing down %d keys took %" PRId64 " milliseconds\n",
    //     par.num, std::chrono::duration_cast<std::chrono::milliseconds>(ed - st).count());
    // printf("Input Enter to Continue...\n");
    // getchar();

    /**
     * Workloads
     */
    std::cout << "Pure range query: " << pureRangeQuery << std::endl;
    if (pureRangeQuery)
    {
        UniformQueryFactory<Key>::getInstance().setup(par.seed, min, max,
                                                      rangeLen,
                                                      {100.0, 0.0, 0.0});
    }
    else
    {
        UniformQueryFactory<Key>::getInstance().setup(par.seed, min, max,
                                                      smallRangeLen, largeRangLen,
                                                      {par.pointReads, par.smallReads, par.largeReads, par.inserts, 0.0});
    }
    memcpy(insertBufferCMD.data(), &insertOP, sizeof(OP));
    std::atomic<bool> stop = false;
    std::thread([&ssl, &opMutex, &stop, par]()
                {
        OP searchOP = OP::SEARCH_OP;
        std::vector<uint8_t> searchBufferCMD;
        std::vector<Datum<Key>> buckets;
        // memcpy(searchBufferCMD.data(), &searchOP, sizeof(OP));
        auto st = std::chrono::high_resolution_clock::now();
        auto ed = st;
        size_t bytesToRead = 0;
        size_t ret = 0;
        while (!stop)
        {
            const std::lock_guard<std::mutex> lock(opMutex);
            SamplePool<Key>::getInstance().sample(searchBufferCMD);
            if (!searchBufferCMD.empty())
            {
                st = std::chrono::high_resolution_clock::now();
                writeToServer(ssl, (uint8_t*)&searchOP, sizeof(searchOP));
                uploadToServer(ssl, searchBufferCMD);
                ret = readFromServer(ssl, (uint8_t*)&bytesToRead, sizeof(bytesToRead));
                if (bytesToRead / sizeof(Datum<Key>) > buckets.size())
                {
                    buckets.resize(bytesToRead / sizeof(Datum<Key>));
                }
                ret = readFromServer(ssl, (uint8_t*)buckets.data(), bytesToRead);
                searchBufferCMD.clear();
                ed = std::chrono::high_resolution_clock::now();
                printf("==========batch delay = %" PRId64 " microseconds========\n",
                    std::chrono::duration_cast<std::chrono::microseconds>(ed - st).count()
                );
            }
            // ed = std::chrono::high_resolution_clock::now();
            // auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(ed - st).count();
            // if (ms < 10)
            // {
            //     std::this_thread::sleep_for(std::chrono::milliseconds(10 - ms));
            // }
        } })
        .detach();

    uint64_t totalTime = 0;

    for (uint32_t i = 0; i < par.queryNum; ++i)
    {
        std::unique_ptr<Query<Key>> curQuery = UniformQueryFactory<Key>::getInstance().createQuery();
        switch (curQuery->getOp())
        {
        case OP::INSERT_OP:
        case OP::DELETE_OP:
            insertBuffer[bufferCnt++].set(curQuery->getQueryData().first, curQuery->getOp(), ++cnt);
            break;
        case OP::SEARCH_OP:
        {
            {
                const std::lock_guard<std::mutex> lock(opMutex);
                SamplePool<Key>::getInstance().insertQuery(
                    curQuery->getQueryData().first,
                    curQuery->getQueryData().second);
            }
            std::this_thread::sleep_for(std::chrono::nanoseconds(1));
            break;
        }
        default:
            throw std::runtime_error("unknown operation");
            break;
        }

        if (bufferCnt == par.bucketCapacity)
        {
            auto st1 = std::chrono::high_resolution_clock::now();
            auto ed1 = st1;
            {
                const std::lock_guard<std::mutex> lock(opMutex);
                st1 = std::chrono::high_resolution_clock::now();
                writeToServer(ssl, (uint8_t *)&insertOP, sizeof(insertOP));
                bufferCnt = 0;
                uploadToServer(ssl, insertBufferCMD);
                updateSamplePool<Key>(ssl);
                readFromServer(ssl, (uint8_t *)&ret, sizeof(ret));
                if (ret != 0)
                {
                    printf(CLIENT_PROXY "Error in setup phase!\n");
                    return ret;
                }
                ed1 = std::chrono::high_resolution_clock::now();
                totalTime += std::chrono::duration_cast<std::chrono::milliseconds>(ed1 - st1).count();
            }
            std::cout << "Inserting took " << std::chrono::duration_cast<std::chrono::milliseconds>(ed1 - st1).count() << " milliseconds" << std::endl;
            // std::cout << "Input q to quit, any other key to continue..." << std::endl;
            // char c = getchar();
            // if (c == 'q')
            // {
            //     exit(0);
            // }
        }
    }
    std::cout << "Total time: " << totalTime << " milliseconds" << std::endl;
    std::this_thread::sleep_for(std::chrono::seconds(1));
    stop = true;
    ret = 0;
done:
    return ret;
}

int createSSLContext(SSL_CTX *&ctx)
{
    if ((ctx = SSL_CTX_new(SSLv23_client_method())) == nullptr)
    {
        printf(CLIENT_PROXY "TLS client: unable to create a new SSL context\n");
        return -1;
    }

    // choose TLSv1.2 by excluding SSLv2, SSLv3 ,TLS 1.0 and TLS 1.1
    SSL_CTX_set_options(ctx, SSL_OP_NO_SSLv2);
    SSL_CTX_set_options(ctx, SSL_OP_NO_SSLv3);
    SSL_CTX_set_options(ctx, SSL_OP_NO_TLSv1);
    SSL_CTX_set_options(ctx, SSL_OP_NO_TLSv1_1);
    // specify the verify_callback for custom verification
    SSL_CTX_set_verify(ctx, SSL_VERIFY_PEER, &verify_callback);
    return 0;
}

int initSSL(SSL *&ssl, SSL_CTX *&ctx, int serversocket)
{
    if ((ssl = SSL_new(ctx)) == nullptr)
    {
        printf(CLIENT_PROXY
               "Unable to create a new SSL connection state object\n");
        return -1;
    }

    // setup ssl socket and initiate TLS connection with TLS server
    SSL_set_fd(ssl, serversocket);
    int error = 0;

    if ((error = SSL_connect(ssl)) != 1)
    {
        printf(
            CLIENT_PROXY "Error: Could not establish an SSL session ret2=%d "
                         "SSL_get_error()=%d\n",
            error,
            SSL_get_error(ssl, error));
        return -1;
    }
    printf(
        CLIENT_PROXY "successfully established TLS channel:%s\n",
        SSL_get_version(ssl));
    return 0;
}

int main(int argc, char **argv)
{
    int ret = 1;
    SSL_CTX *ctx = nullptr;
    SSL *ssl = nullptr;
    ClientParameter par;
    int error = 0;
    int serversocketUpdate = -1;
    int serversocketSearch = -1;
    printf("\nStarting" CLIENT_PROXY "\n\n\n");
    if ((error = parse_arguments(argc, argv, par)) != 0)
    {
        printf(
            CLIENT_PROXY "TLS client:parse input parmeter failed (%d)!\n", error);
        goto done;
    }

    // initialize openssl library and register algorithms
    OpenSSL_add_all_algorithms();
    ERR_load_BIO_strings();
    ERR_load_crypto_strings();
    SSL_load_error_strings();

    if (SSL_library_init() < 0)
    {
        printf(CLIENT_PROXY
               "TLS client: could not initialize the OpenSSL library !\n");
        goto done;
    }

    if (createSSLContext(ctx) != 0)
        goto done;

    serversocketUpdate = create_socket(par.serverName, par.serverPort);
    if (serversocketUpdate == -1)
    {
        printf(
            CLIENT_PROXY
            "create a socket and initate a TCP connect to server: %s:%s "
            "(errno=%d)\n",
            par.serverName,
            par.serverPort,
            errno);
        return -1;
    }
    serversocketSearch = create_socket(par.serverName, par.serverPort);
    if (serversocketSearch == -1)
    {
        printf(
            CLIENT_PROXY
            "create a socket and initate a TCP connect to server: %s:%s "
            "(errno=%d)\n",
            par.serverName,
            par.serverPort,
            errno);
        return -1;
    }
    printf("update socket: %d, search socket: %d\n", serversocketUpdate, serversocketSearch);

    if (initSSL(ssl, ctx, serversocketUpdate) != 0)
        goto done;
    // start the client server communication
    if ((error = sendParameters(ssl, par)) != 0)
    {
        printf(CLIENT_PROXY "Failed: sendParameters (ret=%d)\n", error);
        goto done;
    }

    if ((error = generateQueries<int>(ssl, par)) != 0)
    {
        printf(CLIENT_PROXY "Failed: generateQueries (ret=%d)\n", error);
        goto done;
    }

    // if ((error = generateQueries<Record>(ssl, par)) != 0)
    // {
    //     printf(CLIENT_PROXY "Failed: generateQueries (ret=%d)\n", error);
    //     goto done;
    // }

    // Free the structures we don't need anymore
    ret = 0;
done:
    if (serversocketSearch != -1)
        close(serversocketSearch);
    if (serversocketUpdate != -1)
        close(serversocketUpdate);

    if (ssl)
        SSL_free(ssl);

    if (ctx)
        SSL_CTX_free(ctx);

    printf(CLIENT_PROXY " %s\n", (ret == 0) ? "success" : "failed");
    return (ret);
}
