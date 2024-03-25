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

#include <openssl/evp.h>
#include <openssl/ssl.h>
#include <stdlib.h>
#include <string.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <unistd.h>
#include "common/defs.h"
#include "common/openssl_utility.h"
#include "enclave/DOMerger.hpp"

extern "C"
{
    int set_up_server(void* par_);
    sgx_status_t ocall_close(int* ret, int fd);
    void ocall_getBuckets(uint8_t* bucketIDs_, size_t size_,
        uint8_t** encBuckets, uint32_t encBucketSize, uint32_t* bucketCnt);
};

static sgx_aes_ctr_128bit_key_t* table_key = NULL;

int verify_callback(int preverify_ok, X509_STORE_CTX* ctx);

int create_listener_socket(int port, int& server_socket)
{
    int ret = -1;
    const int reuse = 1;
    struct sockaddr_in addr;
    addr.sin_family = AF_INET;
    addr.sin_port = htons(port);
    addr.sin_addr.s_addr = htonl(INADDR_ANY);

    server_socket = socket(AF_INET, SOCK_STREAM, 0);
    if (server_socket < 0)
    {
        t_print(SERVER_PROXY "socket creation failed\n");
        goto exit;
    }

    if (setsockopt(
        server_socket,
        SOL_SOCKET,
        SO_REUSEADDR,
        (const void*)&reuse,
        sizeof(reuse)) < 0)
    {
        t_print(SERVER_PROXY "setsocket failed \n");
        goto exit;
    }

    if (bind(server_socket, (struct sockaddr*)&addr, sizeof(addr)) < 0)
    {
        t_print(SERVER_PROXY "Unable to bind socket to the port\n");
        goto exit;
    }

    if (listen(server_socket, 20) < 0)
    {
        t_print(SERVER_PROXY "Unable to open socket for listening\n");
        goto exit;
    }
    ret = 0;
exit:
    return ret;
}

template<typename Key>
int handle_queries(
    int& server_socket_fd,
    int& client_socket_fd,
    SSL_CTX*& ssl_server_ctx,
    SSL*& ssl_session,
    std::shared_ptr<BaseIndex>& index)
{
    const size_t bucketCapacity = index->getBucketCapacity();
    const size_t datumSize = index->getDatumSize();
    const size_t bucketSize = bucketCapacity * datumSize;
    const size_t encBucketSize = index->getEncBucketSize();
    int ret = -1;
    t_print(SERVER_PROXY "Bucket size = %d\n", index->getBucketSize());
    size_t bytesToRead = 0;
    std::vector<uint8_t> buffer;
    OP curOP;
    while (1)
    {
        // read operation first 
        if (read_from_session_peer(
            ssl_session, &curOP, sizeof(curOP)) != sizeof(curOP))
        {
            t_print(SERVER_PROXY " Read operation type from client failed\n");
            return ret;
        }

        int ret = 0;
        switch (curOP)
        {
        case OP::INSERT_OP:
        case OP::DELETE_OP:
        {
#ifndef OBLIVIOUSMERGE
            bool sorted = true;
            Key prevKey = std::numeric_limits<Key>::min();
#endif
            if (read_from_session_peer(
                ssl_session, &bytesToRead, sizeof(bytesToRead)) != sizeof(bytesToRead))
            {
                t_print(SERVER_PROXY " Read bytes to read in insertion from client failed\n");
                return ret;
            }

            if (bytesToRead % datumSize != 0)
            {
                t_print(SERVER_PROXY
                    " inconsistent bytesToRead\n");
                return -1;
            }
            if (buffer.size() < bucketSize)
                buffer.resize(bucketSize);

            uint8_t* encBuckets = NULL;
            ocall_resizeMemoryPool(&encBuckets, 1u, bytesToRead);
            uint8_t ctr[16] = { 0 };
            uint8_t* curEncBucket = encBuckets;
            uint32_t dataCnt = bytesToRead / datumSize;
            while (bytesToRead > 0)
            {
                const size_t currentBatchSize = std::min(bucketSize, bytesToRead);
                const uint32_t currentCnt = currentBatchSize / datumSize;
                if (read_from_session_peer(
                    ssl_session, buffer.data(), currentBatchSize) != currentBatchSize)
                {
                    t_print(SERVER_PROXY " Read insertion data from client failed\n");
                    return ret;
                }
                // t_print("Debug remaining bytes = %d, current batch size = %d, sorted = %d\n", bytesToRead, currentBatchSize, sorted);

#ifndef OBLIVIOUSMERGE
                if (sorted)
                {
                    Datum<Key>* curDatum = (Datum<Key>*)(buffer.data());
                    for (uint32_t i = 0; i < currentCnt; i++)
                    {
                        Key curKey = curDatum->getDat();
                        if (curKey < prevKey)
                        {
                            sorted = false;
                            break;
                        }
                        curDatum++;
                        prevKey = curKey;
                    }
                }
#endif

                sgx_aes_ctr_encrypt(table_key,
                    buffer.data(), currentBatchSize,
                    ctr, 128,
                    curEncBucket);
                bytesToRead -= currentBatchSize;
                curEncBucket += currentBatchSize;
            }

#ifndef OBLIVIOUSMERGE
            if (!sorted)
            {
                index->osort(encBuckets, dataCnt);
            }
#endif

            ret = index->insert(encBuckets, dataCnt / bucketCapacity);
            if (ret != 0)
            {
                t_print("Index insert failed!\n");
                return ret;
            }
            break;
        }
        case OP::SEARCH_OP:
        {
            if (read_from_session_peer(
                ssl_session, &bytesToRead, sizeof(bytesToRead)) != sizeof(bytesToRead))
            {
                t_print(SERVER_PROXY " Read bytes to read in search from client failed\n");
                return ret;
            }
            if (bytesToRead > buffer.size())
                buffer.resize(bytesToRead);
            if (read_from_session_peer(
                ssl_session, buffer.data(), bytesToRead) != bytesToRead)
            {
                t_print(SERVER_PROXY " Read search buffer from client failed\n");
                return ret;
            }

            uint8_t* encBuckets = NULL;
            uint32_t bucketCnt = 0;
            ocall_getBuckets(buffer.data(),
                bytesToRead, &encBuckets, encBucketSize, &bucketCnt);

            size_t totalBucketSize = (size_t)encBucketSize * bucketCnt;
            // size_t totalBucketSize = (size_t)bucketSize * bucketCnt;
            // uint8_t* curEncBucket = encBuckets;
            // uint8_t* curBucket = encBuckets;
            // for (uint32_t i = 0; i < bucketCnt; ++i)
            // {
            //     sgx_rijndael128GCM_decrypt(
            //         table_key,
            //         curEncBucket + SGX_AESGCM_IV_SIZE,
            //         bucketSize,
            //         curBucket,
            //         curEncBucket, SGX_AESGCM_IV_SIZE,
            //         NULL, 0,
            //         (sgx_aes_gcm_128bit_tag_t*)(curEncBucket + bucketSize + SGX_AESGCM_IV_SIZE));
            //     curEncBucket += encBucketSize;
            //     curBucket += bucketSize;
            // }
            if (write_to_session_peer(
                ssl_session, &totalBucketSize, sizeof(totalBucketSize)) != 0)
            {
                // ocall_delete(encBuckets);
                t_print(SERVER_PROXY " Write to client failed\n");
                return ret;
            }
            if (write_to_session_peer(
                ssl_session, encBuckets, totalBucketSize) != 0)
            {
                // ocall_delete(encBuckets);
                t_print(SERVER_PROXY " Write to client failed\n");
                return ret;
            }
            // ocall_delete(encBuckets);
            break;
        }
        // case OP::SETUP_OP:
        // {
        //     ret = index->setup(buffer.data() + CMD_LEN);
        //     if (ret != 0)
        //     {
        //         return ret;
        //     }
        //     break;
        // }
        // case OP::PUSH_DOWN_OP:
        // {
        //     ret = index->pushDown();
        //     if (ret != 0)
        //     {
        //         t_print("Index push down failed!\n");
        //         return ret;
        //     }
        //     break;
        // }
        default:
            t_print(SERVER_PROXY " Finished!\n");
            return 0;
        }
    }
    return 0;
}

int handle_parameters(
    int& server_socket_fd,
    int& client_socket_fd,
    SSL_CTX*& ssl_server_ctx,
    SSL*& ssl_session,
    ClientParameter* cPar,
    std::shared_ptr<BaseIndex>& index)
{

    int ret = -1;
    int test_error = 1;

    struct sockaddr_in addr;
    uint len = sizeof(addr);

    client_socket_fd = accept(server_socket_fd, (struct sockaddr*)&addr, &len);

    if (client_socket_fd < 0)
    {
        t_print(SERVER_PROXY "Unable to accept the client request\n");
        goto exitParameter;
    }

    if ((ssl_session = SSL_new(ssl_server_ctx)) == nullptr)
    {
        t_print(SERVER_PROXY
            "Unable to create a new SSL connection state object\n");
        goto exitParameter;
    }

    SSL_set_fd(ssl_session, client_socket_fd);

    test_error = SSL_accept(ssl_session);
    if (test_error <= 0)
    {
        t_print(SERVER_PROXY "SSL handshake failed, error(%d)(%d)\n",
            test_error, SSL_get_error(ssl_session, test_error));
        goto exitParameter;
    }

    if (read_from_session_peer(
        ssl_session, cPar, sizeof(*cPar)) != sizeof(*cPar))
    {
        t_print(SERVER_PROXY "Read from client failed\n");
        goto exitParameter;
    }

    // setup the seed parameter
    ocall_setSeed(cPar->seed);
    genNonce(NULL, 0, cPar->seed); // initialize seed 
    makeBucketID(0, 0, 0, cPar->seed);

    // initialize the index
    switch (cPar->type)
    {
    case Types::Int:
        index = std::make_shared<Index<int>>(table_key, ssl_session, cPar->lambda,
            cPar->epsilon, cPar->delta,
            cPar->bucketCapacity, cPar->storageOverhead,
            cPar->threadNum, cPar->seed, cPar->k);
        // handle communication
        ret = handle_queries<int>(
            server_socket_fd,
            client_socket_fd,
            ssl_server_ctx,
            ssl_session,
            index);
        if (ret != 0)
        {
            t_print(SERVER_PROXY "server communication error %d\n", ret);
            goto exitParameter;
        }
        break;
    case Types::LongLong:
        index = std::make_shared<Index<long long>>(table_key, ssl_session, cPar->lambda,
            cPar->epsilon, cPar->delta,
            cPar->bucketCapacity, cPar->storageOverhead,
            cPar->threadNum, cPar->seed, cPar->k);
        // handle communication
        ret = handle_queries<long long>(
            server_socket_fd,
            client_socket_fd,
            ssl_server_ctx,
            ssl_session,
            index);
        if (ret != 0)
        {
            t_print(SERVER_PROXY "server communication error %d\n", ret);
            goto exitParameter;
        }
        break;
    case Types::Float:
        index = std::make_shared<Index<float>>(table_key, ssl_session, cPar->lambda,
            cPar->epsilon, cPar->delta,
            cPar->bucketCapacity, cPar->storageOverhead,
            cPar->threadNum, cPar->seed, cPar->k);
        ret = handle_queries<float>(
            server_socket_fd,
            client_socket_fd,
            ssl_server_ctx,
            ssl_session,
            index);
        if (ret != 0)
        {
            t_print(SERVER_PROXY "server communication error %d\n", ret);
            goto exitParameter;
        }
        break;
    case Types::Double:
        index = std::make_shared<Index<double>>(table_key, ssl_session, cPar->lambda,
            cPar->epsilon, cPar->delta,
            cPar->bucketCapacity, cPar->storageOverhead,
            cPar->threadNum, cPar->seed, cPar->k);
        ret = handle_queries<double>(
            server_socket_fd,
            client_socket_fd,
            ssl_server_ctx,
            ssl_session,
            index);
        if (ret != 0)
        {
            t_print(SERVER_PROXY "server communication error %d\n", ret);
            goto exitParameter;
        }
        break;
    default:
        throw std::runtime_error("Invalid type");
        break;
    }

    ret = 0;
exitParameter:
    return ret;
}

int set_up_server(void* par_)
{
    int ret = 0;

    ServerParameter* par = (ServerParameter*)par_;
    ClientParameter cPar;
    std::shared_ptr<BaseIndex> index(nullptr);

    int server_socket_fd;
    int client_socket_fd = -1;

    X509* certificate = nullptr;
    EVP_PKEY* pkey = nullptr;
    SSL_CONF_CTX* ssl_confctx = SSL_CONF_CTX_new();

    SSL_CTX* ssl_server_ctx = nullptr;
    SSL* ssl_session = nullptr;
    if ((ssl_server_ctx = SSL_CTX_new(TLS_server_method())) == nullptr)
    {
        t_print(SERVER_PROXY "unable to create a new SSL context\n");
        goto exit;
    }

    if (initalize_ssl_context(ssl_confctx, ssl_server_ctx) != SGX_SUCCESS)
    {
        t_print(SERVER_PROXY "unable to create a initialize SSL context\n ");
        goto exit;
    }
    SSL_CTX_set_verify(ssl_server_ctx, SSL_VERIFY_PEER, &verify_callback);

    if (load_tls_certificates_and_keys(ssl_server_ctx, certificate, pkey) != 0)
    {
        t_print(SERVER_PROXY
            " unable to load certificate and private key on the server\n ");
        goto exit;
    }

    if (create_listener_socket(par->inPort, server_socket_fd) != 0)
    {
        t_print(SERVER_PROXY " unable to create listener socket on the server\n ");
        goto exit;
    }

    // Generate the table key
    table_key = new sgx_aes_ctr_128bit_key_t[SGX_AESCTR_KEY_SIZE];
    sgx_read_rand((uint8_t*)table_key, SGX_AESCTR_KEY_SIZE);

    // Get the client parameters
    ret = handle_parameters(
        server_socket_fd,
        client_socket_fd,
        ssl_server_ctx,
        ssl_session,
        &cPar,
        index);
    if (ret != 0)
    {
        t_print(SERVER_PROXY "server communication error %d\n", ret);
        goto exit;
    }


exit:
    ocall_close(&ret, client_socket_fd); // close the socket connections
    if (ret != 0)
        t_print(SERVER_PROXY "OCALL: error closing client socket\n");
    ocall_close(&ret, server_socket_fd);
    if (ret != 0)
        t_print(SERVER_PROXY "OCALL: error closing server socket\n");

    if (ssl_session)
    {
        SSL_shutdown(ssl_session);
        SSL_free(ssl_session);
    }
    if (ssl_server_ctx)
        SSL_CTX_free(ssl_server_ctx);
    if (ssl_confctx)
        SSL_CONF_CTX_free(ssl_confctx);
    if (certificate)
        X509_free(certificate);
    if (pkey)
        EVP_PKEY_free(pkey);
    return (ret);
}
